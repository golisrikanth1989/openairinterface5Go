/*
 * SPDX-License-Identifier: LicenseRef-CSSL-1.0
 */

/*!
 * \brief Routines for initial UE synchronization procedure (PSS,SSS,PBCH and frame format detection)
 */
#include "PHY/defs_nr_UE.h"
#include "PHY/MODULATION/modulation_UE.h"
#include "nr_transport_proto_ue.h"
#include "PHY/NR_UE_ESTIMATION/nr_estimation.h"
#include "SCHED_NR_UE/defs.h"
#include "common/utils/nr/nr_common.h"

#include <math.h>

#include "PHY/NR_REFSIG/pss_nr.h"
#include "PHY/NR_REFSIG/sss_nr.h"
#include "PHY/NR_REFSIG/nr_refsig.h"
#include "PHY/MODULATION/nr_modulation.h"
#include "PHY/TOOLS/tools_defs.h"
#include "nr-uesoftmodem.h"
#include "nfapi/open-nFAPI/nfapi/public_inc/fapi_nr_ue_interface.h"

//#define DEBUG_INITIAL_SYNCH
#define DUMP_PBCH_CH_ESTIMATES 0

// structure used for multiple SSB detection
typedef struct NR_UE_SSB {
  uint i_ssb; // i_ssb between 0 and 7 (it corresponds to ssb_index only for Lmax=4,8)
  uint n_hf; // n_hf = 0,1 for Lmax =4 or n_hf = 0 for Lmax =8,64
  double metric; // metric to order SSB hypothesis
} NR_UE_SSB;

static int ssb_sort(const void *a, const void *b)
{
  return ((NR_UE_SSB *)b)->metric - ((NR_UE_SSB *)a)->metric;
}

static bool nr_pbch_detection(const UE_nr_rxtx_proc_t *proc,
                              const NR_DL_FRAME_PARMS *frame_parms,
                              int Nid_cell,
                              int pbch_initial_symbol,
                              int ssb_start_subcarrier,
                              int *half_frame_bit,
                              int *ssb_index,
                              int *symbol_offset,
                              fapiPbch_t *result,
                              const c16_t rxdataF[NR_N_SYMBOLS_SSB][frame_parms->nb_antennas_rx][frame_parms->ofdm_symbol_size])
{
  const int N_L = (frame_parms->Lmax == 4) ? 4 : 8;
  const int N_hf = (frame_parms->Lmax == 4) ? 2 : 1;
  NR_UE_SSB best_ssb[N_L * N_hf];
  NR_UE_SSB *current_ssb = best_ssb;
  // loops over possible pbch dmrs cases to retrieve best estimated i_ssb (and n_hf for Lmax=4) for multiple ssb detection
  for (int hf = 0; hf < N_hf; hf++) {
    for (int l = 0; l < N_L; l++) {
      // computing correlation between received DMRS symbols and transmitted sequence for current i_ssb and n_hf
      cd_t cumul = {0};
      for (int i = pbch_initial_symbol; i < pbch_initial_symbol + 3; i++) {
        c32_t meas = nr_pbch_dmrs_correlation(frame_parms,
                                              i,
                                              i - pbch_initial_symbol,
                                              Nid_cell,
                                              ssb_start_subcarrier,
                                              nr_gold_pbch(frame_parms->Lmax, Nid_cell, hf, l),
                                              rxdataF[i]);
        csum(cumul, cumul, meas);
      }
      *current_ssb = (NR_UE_SSB){.i_ssb = l, .n_hf = hf, .metric = squaredMod(cumul)};
      current_ssb++;
    }
  }
  qsort(best_ssb, N_L * N_hf, sizeof(NR_UE_SSB), ssb_sort);

  const int nb_ant = frame_parms->nb_antennas_rx;
  const int estimateSz = frame_parms->ofdm_symbol_size;
  for (NR_UE_SSB *ssb = best_ssb; ssb < best_ssb + N_L * N_hf; ssb++) {
    // computing channel estimation for selected best ssb
    int16_t pbch_e_rx[NR_POLAR_PBCH_E];

    uint8_t log2_maxh = 0;
    for (int i = pbch_initial_symbol; i < pbch_initial_symbol + 3; i++) {
      __attribute__((aligned(32))) c16_t dl_ch_estimates[nb_ant][estimateSz];
      for (int aarx = 0; aarx < nb_ant; aarx++) {
        nr_pbch_channel_estimation(frame_parms,
                                   NULL,
                                   dl_ch_estimates[aarx],
                                   proc,
                                   i - pbch_initial_symbol,
                                   ssb->i_ssb,
                                   ssb->n_hf,
                                   ssb_start_subcarrier,
                                   rxdataF[i][aarx],
                                   false,
                                   Nid_cell);
      }
      if (DUMP_PBCH_CH_ESTIMATES) {
        char varName[30] = "";
        snprintf(varName, sizeof(varName), "pbch_ch_estimates_symbol_%d", i);
        LOG_MM("pbch_ch_estimates", varName, dl_ch_estimates, nb_ant * estimateSz, 1, 1);
      }
      nr_generate_pbch_llr(NULL,
                           proc,
                           frame_parms,
                           i,
                           ssb->i_ssb,
                           Nid_cell,
                           ssb_start_subcarrier,
                           rxdataF[i],
                           dl_ch_estimates,
                           pbch_e_rx,
                           &log2_maxh);
    }

    if (0
        == nr_pbch_decode(NULL,
                          frame_parms,
                          proc,
                          ssb->i_ssb,
                          Nid_cell,
                          pbch_e_rx,
                          half_frame_bit,
                          ssb_index,
                          symbol_offset,
                          result)) {
      LOG_A(PHY, "Initial sync: pbch decoded sucessfully, ssb index %d\n", *ssb_index);
      return true;
    }
  }

  LOG_W(PHY, "Initial sync: pbch not decoded, ssb index %d\n", frame_parms->ssb_index);
  return false;
}

static void compensate_freq_offset(c16_t **x, const int nb_antennas_rx, const int x_len, const int offset, const int sampling_rate)
{
  double s_time = 1.0 / sampling_rate; // sampling time
  double off_angle = -2 * M_PI * s_time * (offset); // offset rotation angle compensation per sample

  for (int n = 0; n < x_len; n++) {
    for (int ar = 0; ar < nb_antennas_rx; ar++) {
      const double re = x[ar][n].r;
      const double im = x[ar][n].i;
      x[ar][n].r = (short)(round(re * cos(n * off_angle) - im * sin(n * off_angle)));
      x[ar][n].i = (short)(round(re * sin(n * off_angle) + im * cos(n * off_angle)));
    }
  }
}

/*
 * Van de Beek / Sandell / Boerjesson cyclic-prefix fractional frequency offset estimator.
 * (J.J. van de Beek et al., "ML estimation of time and frequency offset in OFDM systems",
 *  IEEE Trans. Signal Process., vol. 45, no. 7, 1997.)
 *
 * The CP of each OFDM symbol is a copy of that symbol's last Ncp samples, i.e. the same
 * samples seen N samples apart.  A CFO eps (in units of the SCS) rotates the later copy by
 *      r(n + N) = r(n) * exp(j*2*pi*eps)
 * so the conjugate product has the data magnitude (real, positive) and the CFO phase only:
 *      conj(r(n)) * r(n + N) = |s(n)|^2 * exp(j*2*pi*eps)
 * Accumulating over the Ncp CP samples of the 4 SSB symbols (PSS, PBCH, SSS, PBCH) and
 * taking the angle gives 2*pi*eps, hence the division by 2*pi (OAI's PSS half-symbol
 * estimator uses a lag of N/2, giving pi*eps and a division by pi instead).
 *
 * Lag N -> unambiguous range |eps| < 0.5, i.e. +-SCS/2 (half that of the PSS estimator),
 * but twice the phase per Hz, so roughly twice the accuracy at the same SNR.
 *
 * Returns the estimate in Hz.  rxdata must still carry the offset to be measured; any
 * offset already compensated in rxdata has to be added back by the caller.
 */
static int estimate_ffo_vdb_cp(c16_t **rxdata,
                               int nb_antennas_rx,
                               int ssb_time_offset,
                               int ofdm_symbol_size,
                               int nb_prefix_samples,
                               int subcarrier_spacing)
{
  const int N = ofdm_symbol_size;
  const int Ncp = nb_prefix_samples;
  const int L = N + Ncp; // one OFDM symbol, CP included
  double acc_re = 0.0, acc_im = 0.0;

  for (int ar = 0; ar < nb_antennas_rx; ar++) {
    for (int m = 0; m < NR_N_SYMBOLS_SSB; m++) {
      const c16_t *sym = &rxdata[ar][ssb_time_offset + m * L];
      for (int n = 0; n < Ncp; n++) {
        const c16_t a = sym[n]; // CP sample
        const c16_t b = sym[n + N]; // the same sample N later (end of the useful part)
        acc_re += (double)a.r * b.r + (double)a.i * b.i; // Re{conj(a) * b}
        acc_im += (double)a.r * b.i - (double)a.i * b.r; // Im{conj(a) * b}
      }
    }
  }

  if (acc_re == 0.0 && acc_im == 0.0)
    return 0;

  return (int)lround(subcarrier_spacing * atan2(acc_im, acc_re) / (2 * M_PI));
}

/*
 * Optional CSV log of both FFO estimators, enabled with --ffo-log <file>.
 * One row per PSS candidate per sync attempt; the file is opened on first use and
 * a header row is written when the file is created.
 */
static FILE *ffo_log_fd = NULL;
static pthread_mutex_t ffo_log_lock = PTHREAD_MUTEX_INITIALIZER;

static void ffo_log_row(bool fo_flag,
                        int nid2,
                        int pss_pos,
                        int peak_db,
                        int avg_db,
                        int ffo_oai_hz,
                        int ffo_vdb_hz,
                        int ffo_pre_hz,
                        int ffo_used_hz,
                        int use_vdb,
                        int scs_hz,
                        int fft_size,
                        int ncp)
{
  const char *path = get_nrUE_params()->ffo_log;
  if (path == NULL)
    return;

  // Both estimators are gated on fo_flag, so with it off every estimate column
  // would be 0. Log only rows where a value was actually computed. This also
  // drops the neighbour-cell search, which calls nr_search_ssb_common() with
  // fo_flag = false (nr_ue_measurements.c) and runs continuously once the UE is
  // connected -- it would otherwise bury the real sync rows under hundreds of
  // empty ones.
  if (!fo_flag)
    return;

  pthread_mutex_lock(&ffo_log_lock);
  if (ffo_log_fd == NULL) {
    ffo_log_fd = fopen(path, "a");
    if (ffo_log_fd == NULL) {
      LOG_E(PHY, "could not open FFO log file %s\n", path);
      pthread_mutex_unlock(&ffo_log_lock);
      return;
    }
    if (ftell(ffo_log_fd) == 0)
      fprintf(ffo_log_fd,
              "t_us,nid2,pss_pos,peak_db,avg_db,corr_db,"
              "ffo_oai_hz,ffo_vdb_hz,ffo_pre_hz,ffo_used_hz,"
              "eps_oai,eps_vdb,method,scs_hz,fft_size,ncp\n");
  }

  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  const uint64_t t_us = (uint64_t)ts.tv_sec * 1000000 + ts.tv_nsec / 1000;

  fprintf(ffo_log_fd,
          "%lu,%d,%d,%d,%d,%d,%d,%d,%d,%d,%.6f,%.6f,%s,%d,%d,%d\n",
          (unsigned long)t_us,
          nid2,
          pss_pos,
          peak_db,
          avg_db,
          peak_db - avg_db,
          ffo_oai_hz,
          ffo_vdb_hz,
          ffo_pre_hz,
          ffo_used_hz,
          scs_hz ? (double)ffo_oai_hz / scs_hz : 0.0,
          scs_hz ? (double)ffo_vdb_hz / scs_hz : 0.0,
          use_vdb ? "vdb_cp" : "oai_pss",
          scs_hz,
          fft_size,
          ncp);
  fflush(ffo_log_fd); // flush per row: the UE may be killed with Ctrl-C mid-run
  pthread_mutex_unlock(&ffo_log_lock);
}

/* rxdataF should be 16 bytes aligned */
static void generate_table(nr_ssb_search_params_t *params,
                           c16_t timeshift_symbol_rotation[params->ofdm_symbol_size],
                           c16_t symbol_rotation[224])
{
  init_timeshift_rotation(params->ofdm_symbol_size,
                          params->N_RB_DL * NR_NB_SC_PER_RB,
                          params->nb_prefix_samples,
                          params->ofdm_offset_divisor,
                          timeshift_symbol_rotation);
  perform_symbol_rotation(params->symbols_per_slot * params->slots_per_frame / 10,
                          params->numerology_index,
                          params->dl_CarrierFreq,
                          symbol_rotation);
}

static void do_time_to_freq(nr_ssb_search_params_t *params, uint32_t sample_offset)
{
  c16_t timeshift_symbol_rotation[params->ofdm_symbol_size];
  c16_t symbol_rotation[224];
  generate_table(params, timeshift_symbol_rotation, symbol_rotation);

  c16_t(*rxdataF)[params->nb_antennas_rx][params->ofdm_symbol_size] =
      (c16_t(*)[params->nb_antennas_rx][params->ofdm_symbol_size])params->rxdataF;
  dft_size_idx_t dftsize = get_dft(params->ofdm_symbol_size);

  for (int symb = 0; symb < NR_N_SYMBOLS_SSB; symb++) {
    // For Sidelink 16 frames worth of samples is processed to find SSB, for 5G-NR 2.
    unsigned int rx_offset = sample_offset + params->nb_prefix_samples;
    rx_offset += symb * (params->nb_prefix_samples + params->ofdm_symbol_size);
    // use OFDM symbol from within 1/8th of the CP to avoid ISI
    rx_offset -= params->nb_prefix_samples / params->ofdm_offset_divisor;
    for (unsigned char aa = 0; aa < params->nb_antennas_rx; aa++) {
      c16_t *rxF = rxdataF[symb][aa];
      // OFDM Demod
      dft(dftsize, (int16_t *)&params->rxdata[aa][rx_offset], (int16_t *)rxF, 1);
      // FFT-shift
      fftshift_inplace(rxF, params->N_RB_DL * NR_NB_SC_PER_RB, params->ofdm_symbol_size);
      // Phase compensation
      apply_nr_rotation_symbol_fftshifted_RX(params->symbols_per_slot,
                                             params->slots_per_subframe,
                                             timeshift_symbol_rotation,
                                             rxF,
                                             symbol_rotation,
                                             params->N_RB_DL,
                                             0,
                                             symb);
    }
  }
}

/*
 * Common SSB search function used by both initial sync and neighbor cell search
 */
bool nr_search_ssb_common(nr_ssb_search_params_t *params)
{
  const uint32_t pssTime_sz = params->ofdm_symbol_size;
  c16_t(*pssTime)[pssTime_sz] = (c16_t(*)[pssTime_sz])params->pssTime;

  // Perform PSS search
  pss_search_t p_pss = (pss_search_t){.rxdata = params->rxdata,
                                      .nb_antennas_rx = params->nb_antennas_rx,
                                      .rxdata_length = params->rxdata_size,
                                      .ofdm_symbol_size = params->ofdm_symbol_size,
                                      .nb_prefix_samples = params->nb_prefix_samples,
                                      .subcarrier_spacing = params->subcarrier_spacing,
                                      .fo_flag = params->fo_flag,
                                      .target_Nid_cell = params->target_nid_cell,
                                      .pssTime = (c16_t *)pssTime};
  nr_pss_info_t pss_info = pss_search_time_nr(&p_pss);

  // This is the frequency offset that will be applied in the compensation,
  // and it takes into account the values already applied previously during the loop.
  int f_off_to_comp = 0;
  // Net frequency shift currently applied to params->rxdata (0 on the first iteration).
  // The CP estimator measures what is left in the samples, so this has to be added back.
  int f_off_already_applied = 0;
  const bool use_vdb = get_nrUE_params()->vbk_cp_ffo != 0;
  const bool ffo_logging = get_nrUE_params()->ffo_log != NULL;

  for (int p = 0; p < NUMBER_PSS_SEQUENCE; p++) {
    pss_detection_result_t *pss_res = &pss_info.pss_elem_info[p];
    if (!pss_res->success)
      continue;

    int freq_offset_pss = pss_res->freq_offset;
    int sync_pos = pss_res->pos;

    const int ssb_time_offset = sync_pos - params->nb_prefix_samples;

#ifdef DEBUG_INITIAL_SYNCH
    LOG_I(PHY, "Initial sync : Estimated PSS position %d, Nid2 %d, ssb time offset %d\n", sync_pos, p, ssb_time_offset);
#endif

    // Check that SSB fits within buffer
    if (ssb_time_offset + NR_N_SYMBOLS_SSB * (params->ofdm_symbol_size + params->nb_prefix_samples) >= params->rxdata_size) {
      LOG_D(PHY,
            "SSB extends beyond buffer boundary (sync_pos %d, ssb_time_offset %d, buffer_size %d)\n",
            sync_pos,
            ssb_time_offset,
            params->rxdata_size);
      return false;
    }

    // Fractional frequency offset: OAI's PSS half-symbol estimator (default) or the
    // Van de Beek cyclic-prefix estimator (--vbk-cp-ffo). Both are computed whenever
    // --ffo-log is given, so that a single run yields directly comparable data.
    int freq_offset_vdb = 0;
    if (params->fo_flag && (use_vdb || ffo_logging)) {
      // the SSB is known to fit in the buffer at this point (checked just above)
      freq_offset_vdb = f_off_already_applied
                        + estimate_ffo_vdb_cp(params->rxdata,
                                              params->nb_antennas_rx,
                                              ssb_time_offset,
                                              params->ofdm_symbol_size,
                                              params->nb_prefix_samples,
                                              params->subcarrier_spacing);
    }
    const int freq_offset_used = use_vdb ? freq_offset_vdb : freq_offset_pss;

    ffo_log_row(params->fo_flag,
                pss_res->nid2,
                sync_pos,
                pss_res->peak,
                pss_res->avg,
                freq_offset_pss,
                freq_offset_vdb,
                f_off_already_applied,
                freq_offset_used,
                use_vdb,
                params->subcarrier_spacing,
                params->ofdm_symbol_size,
                params->nb_prefix_samples);

    LOG_D(PHY,
          "FFO estimates: OAI PSS half-symbol %d Hz, Van de Beek CP %d Hz, using %s -> %d Hz\n",
          freq_offset_pss,
          freq_offset_vdb,
          use_vdb ? "vdb_cp" : "oai_pss",
          freq_offset_used);

    // the rest of the chain (SSS refinement, accumulated freqOffset) uses pss_res->freq_offset
    pss_res->freq_offset = freq_offset_used;

    // Apply frequency offset compensation if requested
    if (params->apply_freq_offset && freq_offset_used != 0) {
      f_off_to_comp += freq_offset_used;
      compensate_freq_offset(params->rxdata, params->nb_antennas_rx, params->rxdata_size, f_off_to_comp, params->sampling_rate);
      f_off_to_comp *= -1;
      f_off_already_applied = freq_offset_used;
    }

    // Extract SSB symbols to frequency domain
    // Symbol ordering: 0=PSS, 1=PBCH, 2=SSS, 3=PBCH
    do_time_to_freq(params, ssb_time_offset);

    // Perform SSS detection
    nr_sss_params_t p_sss = (nr_sss_params_t){.nb_antennas_rx = params->nb_antennas_rx,
                                              .samples_per_slot_wCP = params->samples_per_slot_wCP,
                                              .ofdm_symbol_size = params->ofdm_symbol_size,
                                              .ssb_start_subcarrier = params->ssb_start_subcarrier,
                                              .subcarrier_spacing = params->subcarrier_spacing,
                                              .exclude_nid_cells = params->exclude_nid_cells,
                                              .num_exclude_nid_cells = params->num_exclude_nid_cells};

    c16_t(*rxdataF)[params->nb_antennas_rx][params->ofdm_symbol_size] =
        (c16_t(*)[params->nb_antennas_rx][params->ofdm_symbol_size])params->rxdataF;
    params->sss_res = rx_sss_nr(&p_sss, pss_res, -1, rxdataF);

    if (!params->sss_res.success || params->sss_res.nid_cell < 0) {
      continue;
    }

    params->pss_res = *pss_res;
    return true;
  }

  return false;
}

static void nr_scan_ssb(void *arg)
{
  /*   Initial synchronisation
   *
   *                                 1 radio frame = 10 ms
   *     <--------------------------------------------------------------------------->
   *     -----------------------------------------------------------------------------
   *     |                                 Received UE data buffer                    |
   *     ----------------------------------------------------------------------------
   *                     --------------------------
   *     <-------------->| pss | pbch | sss | pbch |
   *                     --------------------------
   *          sync_pos            SS/PBCH block
   */

  nr_ue_ssb_scan_t *ssbInfo = (nr_ue_ssb_scan_t *)arg;
  c16_t **rxdata = ssbInfo->rxdata;
  const NR_DL_FRAME_PARMS *fp = ssbInfo->fp;

  // Generate PSS time signal for this GSCN.
  __attribute__((aligned(32))) c16_t pssTime[NUMBER_PSS_SEQUENCE][fp->ofdm_symbol_size];
  const int pss_sequence = get_softmodem_params()->sl_mode == 0 ? NUMBER_PSS_SEQUENCE : NUMBER_PSS_SEQUENCE_SL;
  for (int nid2 = 0; nid2 < pss_sequence; nid2++)
    generate_pss_nr_time(fp->ofdm_symbol_size, fp->first_carrier_offset, nid2, ssbInfo->gscnInfo.ssbFirstSC, pssTime[nid2]);

  __attribute__((aligned(32))) c16_t rxdataF[NR_N_SYMBOLS_SSB][fp->nb_antennas_rx][fp->ofdm_symbol_size];

  // initial sync performed on two successive frames, if pbch passes on first frame, no need to process second frame
  // only one frame is used for simulation tools
  if (ssbInfo->freqOffset)
    compensate_freq_offset(rxdata, fp->nb_antennas_rx, ssbInfo->rxdata_sz, ssbInfo->freqOffset, fp->samples_per_subframe * 1000);

  nr_ssb_search_params_t search_params = {
      .dl_CarrierFreq = fp->dl_CarrierFreq,
      .sampling_rate = fp->samples_per_subframe * 1000,
      .slots_per_frame = fp->slots_per_frame,
      .slots_per_subframe = fp->slots_per_subframe,
      .numerology_index = fp->numerology_index,
      .ofdm_symbol_size = fp->ofdm_symbol_size,
      .ofdm_offset_divisor = fp->ofdm_offset_divisor,
      .nb_antennas_rx = fp->nb_antennas_rx,
      .symbols_per_slot = fp->symbols_per_slot,
      .N_RB_DL = fp->N_RB_DL,
      .rxdata_size = ssbInfo->rxdata_sz,
      .rxdata = rxdata,
      .nb_prefix_samples = fp->nb_prefix_samples,
      .nb_prefix_samples0 = fp->nb_prefix_samples0,
      .ssb_start_subcarrier = ssbInfo->gscnInfo.ssbFirstSC,
      .subcarrier_spacing = fp->subcarrier_spacing,
      .samples_per_slot_wCP = fp->samples_per_slot_wCP,
      .target_nid_cell = ssbInfo->targetNidCell,
      .exclude_nid_cells = NULL, // No exclusion for initial sync
      .num_exclude_nid_cells = 0,
      .apply_freq_offset = ssbInfo->foFlag,
      .fo_flag = ssbInfo->foFlag,
      .rxdataF = rxdataF,
      .pssTime = pssTime,
  };

  ssbInfo->syncRes.cell_detected = nr_search_ssb_common(&search_params);

  ssbInfo->pssCorrAvgPower = search_params.pss_res.avg;
  ssbInfo->pssCorrPeakPower = search_params.pss_res.peak;
  ssbInfo->ssbOffset = search_params.pss_res.pos - search_params.nb_prefix_samples;
  ssbInfo->nidCell = search_params.sss_res.nid_cell;

#ifdef DEBUG_INITIAL_SYNCH
  LOG_I(PHY,
        "TDD Normal prefix: sss detection result; %d, CellId %d metric %d, phase %d, measured offset %d\n",
        ssbInfo->syncRes.cell_detected,
        ssbInfo->nidCell,
        sss_metric,
        sss_phase,
        ssbInfo->syncRes.rx_offset);
#endif
  ssbInfo->freqOffset += search_params.pss_res.freq_offset + search_params.sss_res.freq_offset;

  if (ssbInfo->syncRes.cell_detected) { // we got sss channel
    ssbInfo->syncRes.cell_detected = nr_pbch_detection(ssbInfo->proc,
                                                       ssbInfo->fp,
                                                       ssbInfo->nidCell,
                                                       1,
                                                       ssbInfo->gscnInfo.ssbFirstSC,
                                                       &ssbInfo->halfFrameBit,
                                                       &ssbInfo->ssbIndex,
                                                       &ssbInfo->symbolOffset,
                                                       &ssbInfo->pbchResult,
                                                       rxdataF); // start pbch detection at first symbol after pss
    if (ssbInfo->syncRes.cell_detected) {
      uint32_t rsrp_avg = nr_ue_calculate_ssb_rsrp(ssbInfo->fp, rxdataF[2], ssbInfo->gscnInfo.ssbFirstSC);
      int rsrp_db_per_re = 10 * log10(rsrp_avg);
      ssbInfo->adjust_rxgain = TARGET_RX_POWER - rsrp_db_per_re;
      LOG_I(PHY, "pbch rx ok. rsrp:%d dB/RE, adjust_rxgain:%d dB\n", rsrp_db_per_re, ssbInfo->adjust_rxgain);
    }
  }
  completed_task_ans(ssbInfo->ans);
}

nr_initial_sync_t nr_initial_sync(UE_nr_rxtx_proc_t *proc,
                                  PHY_VARS_NR_UE *ue,
                                  int input_sz,
                                  c16_t **input,
                                  nr_gscn_info_t gscnInfo[MAX_GSCN_BAND],
                                  int numGscn)
{
  NR_DL_FRAME_PARMS *fp = &ue->frame_parms;

  // Perform SSB scanning in parallel. One GSCN per thread.
  LOG_I(NR_PHY,
        "Starting cell search with center freq: %ld, bandwidth: %d. Scanning for %d number of GSCN.\n",
        fp->dl_CarrierFreq,
        fp->N_RB_DL,
        numGscn);
  DevAssert(numGscn);
  task_ans_t ans;
  init_task_ans(&ans, numGscn);
  nr_ue_ssb_scan_t ssb_info[numGscn];
  for (int s = 0; s < numGscn; s++) {
    nr_ue_ssb_scan_t *ssbInfo = &ssb_info[s];
    *ssbInfo = (nr_ue_ssb_scan_t){.gscnInfo = gscnInfo[s],
                                  .fp = &ue->frame_parms,
                                  .proc = proc,
                                  .syncRes.cell_detected = false,
                                  .foFlag = ue->UE_fo_compensation,
                                  .freqOffset = ue->initial_fo,
                                  .targetNidCell = ue->target_Nid_cell};
    ssbInfo->rxdata = malloc16_clear(fp->nb_antennas_rx * sizeof(c16_t *));
    for (int ant = 0; ant < fp->nb_antennas_rx; ant++) {
      ssbInfo->rxdata[ant] = malloc16(sizeof(c16_t) * input_sz);
      memcpy(ssbInfo->rxdata[ant], input[ant], sizeof(c16_t) * input_sz);
      ssbInfo->rxdata_sz = input_sz;
    }
    LOG_I(NR_PHY,
          "Scanning GSCN: %d, with SSB offset: %d, SSB Freq: %lf\n",
          ssbInfo->gscnInfo.gscn,
          ssbInfo->gscnInfo.ssbFirstSC,
          ssbInfo->gscnInfo.ssRef);
    ssbInfo->ans = &ans;
    task_t t = {.func = nr_scan_ssb, .args = ssbInfo};
    pushTpool(&get_nrUE_params()->Tpool, t);
  }

  // Collect the scan results
  nr_ue_ssb_scan_t *res = NULL;
  join_task_ans(&ans);
  for (int i = 0; i < numGscn; i++) {
    nr_ue_ssb_scan_t *ssbInfo = &ssb_info[i];
    if (ssbInfo->syncRes.cell_detected) {
      LOG_I(NR_PHY,
            "Cell Detected with GSCN: %d, SSB SC offset: %d, SSB Ref: %lf, PSS Corr peak: %d dB, PSS Corr Average: %d\n",
            ssbInfo->gscnInfo.gscn,
            ssbInfo->gscnInfo.ssbFirstSC,
            ssbInfo->gscnInfo.ssRef,
            ssbInfo->pssCorrPeakPower,
            ssbInfo->pssCorrAvgPower);
      // take the first cell detected
      if (!res)
        res = ssbInfo;
    }
    for (int ant = 0; ant < fp->nb_antennas_rx; ant++) {
      free(ssbInfo->rxdata[ant]);
    }
    free(ssbInfo->rxdata);
    ssbInfo->rxdata = NULL;
  }

  // Set globals based on detected cell
  if (res) {
    fp->Nid_cell = res->nidCell;
    fp->ssb_start_subcarrier = res->gscnInfo.ssbFirstSC;
    fp->half_frame_bit = res->halfFrameBit;
    fp->ssb_index = res->ssbIndex;
    ue->symbol_offset = res->symbolOffset;
    ue->common_vars.freq_offset = res->freqOffset;
    ue->adjust_rxgain = res->adjust_rxgain;
  }

  // In initial sync, we indicate PBCH to MAC after the scan is complete.
  if (ue->if_inst && ue->if_inst->dl_indication) {
    fapi_nr_rx_indication_t rx_ind;
    rx_ind.number_pdus = 0;
    nr_fill_rx_indication(&rx_ind, FAPI_NR_RX_PDU_TYPE_SSB, ue, 0, 0, NULL, proc, res ? (void *)&res->pbchResult : NULL);
    nr_downlink_indication_t dl_indication = (nr_downlink_indication_t){
        .gNB_index = proc->gNB_id,
        .module_id = ue->Mod_id,
        .cc_id = ue->CC_id,
        .hfn = proc->hfn_rx,
        .frame = proc->frame_rx,
        .slot = proc->nr_slot_rx,
        .rx_ind = &rx_ind,
    };
    ue->if_inst->dl_indication(&dl_indication);
  }

  LOG_D(PHY, "nr_initial sync ue RB_DL %d\n", fp->N_RB_DL);

  if (res) {
    // sync at symbol ue->symbol_offset
    // computing the offset wrt the beginning of the frame
    int mu = fp->numerology_index;
    // number of symbols with different prefix length
    // every 7*(1<<mu) symbols there is a different prefix length (38.211 5.3.1)
    int n_symb_prefix0 = (res->symbolOffset / (7 * (1 << mu))) + 1;
    const int sync_pos_frame = n_symb_prefix0 * (fp->ofdm_symbol_size + fp->nb_prefix_samples0)
                               + (res->symbolOffset - n_symb_prefix0) * (fp->ofdm_symbol_size + fp->nb_prefix_samples);

    res->syncRes.rx_offset = res->ssbOffset - sync_pos_frame;

    LOG_I(PHY, "[UE%d] In synch, rx_offset %d samples\n", ue->Mod_id, res->syncRes.rx_offset);
    LOG_I(PHY, "[UE %d] Measured Carrier Frequency offset %d Hz\n", ue->Mod_id, res->freqOffset);
    LOG_A(PHY, "Initial sync successful, PCI: %d\n", fp->Nid_cell);
    return res->syncRes;
  } else {
#ifdef DEBUG_INITIAL_SYNC
    LOG_I(PHY,"[UE%d] Initial sync : PBCH not ok\n",ue->Mod_id);
    LOG_I(PHY, "[UE%d] Initial sync : Estimated PSS position %d, Nid2 %d\n", ue->Mod_id, sync_pos, ue->common_vars.nid2);
    LOG_I(PHY,"[UE%d] Initial sync : Estimated Nid_cell %d, Frame_type %d\n",ue->Mod_id,
          fp->Nid_cell,fp->frame_type);
    LOG_I(PHY, "[UE%d] Initial sync failed : Estimated power: %d dB\n", ue->Mod_id, ue->measurements.rx_power_avg_dB[0]);
#endif
    // gain control
    // we are not synched, so we cannot use rssi measurement (which is based on channel estimates)
    int rx_power = 0;

    // do a measurement on the best guess of the PSS
    // for (aarx=0; aarx<frame_parms->nb_antennas_rx; aarx++)
    //  rx_power += signal_energy(&ue->common_vars.rxdata[aarx][sync_pos2],
    //			frame_parms->ofdm_symbol_size+frame_parms->nb_prefix_samples);

    /*
    // do a measurement on the full frame
    for (aarx=0; aarx<frame_parms->nb_antennas_rx; aarx++)
    rx_power += signal_energy(&ue->common_vars.rxdata[aarx][0],
    frame_parms->samples_per_subframe*10);
    */

    // we might add a low-pass filter here later
    ue->measurements.rx_power_avg[0] = rx_power / fp->nb_antennas_rx;
    ue->measurements.rx_power_avg_dB[0] = dB_fixed(ue->measurements.rx_power_avg[0]);
    return (nr_initial_sync_t){.cell_detected = false};
  }
}
