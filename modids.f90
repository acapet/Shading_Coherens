!
! Copyright 2013 RBINS-MUMM
!
! Licensed under the EUPL, Version 1.1 or - as soon they will be approved by
! the European Commission - subsequent versions of the EUPL (the "Licence");
! You may not use this work except in compliance with the Licence.
! You may obtain a copy of the Licence at:
!
! http://ec.europa.eu/idabc/eupl
!
! Unless required by the applicable law or agreed to in writing, software
! distributed under the Licence is distributed on an "AS IS" basis,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the Licence for the specific language governing permissions and
! limitations under the Licence.

MODULE modids
!************************************************************************
!
! *modids* Model variable key ids
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)modids.f90  V2.11.3
!
! $Date: 2021-11-28 10:15:48 +0100 (So, 28 Nov 2021) $
!
! $Revision: 1435 $
!
! Description -
!
! Reference -
!
!************************************************************************
!

USE syspars

IMPLICIT NONE

!---model grid
INTEGER, PARAMETER :: &
& iarr_alphatc_fld = 1, iarr_alphatu_fld = 2, iarr_alphatv_fld = 3, &
& iarr_coriolatu = 4, iarr_coriolatv = 5, iarr_delxatc = 6, iarr_delxatu = 7, &
& iarr_delxatuv = 8, iarr_delxatv = 9, iarr_delyatc = 10, iarr_delyatu = 11, &
& iarr_delyatuv = 12, iarr_delyatv = 13, iarr_delzatc = 14, iarr_delzatu = 15, &
& iarr_delzatuv = 16, iarr_delzatuw = 17, iarr_delzatv = 18, &
& iarr_delzatvw = 19, iarr_delzatw = 20, iarr_dryfac = 21, iarr_gaccatc = 22, &
& iarr_gaccatu = 23, iarr_gaccatv = 24, iarr_gangleatc = 25, iarr_garea = 26, &
& iarr_gdelxglb = 27, iarr_gdelyglb = 28, iarr_gscoordatc = 29, &
& iarr_gscoordatu = 30, iarr_gscoordatuvw = 31, iarr_gscoordatuw = 32, &
& iarr_gscoordatv = 33, iarr_gscoordatvw = 34, iarr_gscoordatw = 35, &
& iarr_gscoordglb = 36, iarr_gsigcoordatc = 37, iarr_gsigcoordatw = 38, &
& iarr_gxcoord = 39, iarr_gxcoordglb = 40, iarr_gxlon = 41, iarr_gycoord = 42, &
& iarr_gycoordglb = 43, iarr_gylat = 44, iarr_indexobu = 45, &
& iarr_indexobuprocs = 46, iarr_indexobv = 47, iarr_indexobvprocs = 48, &
& iarr_indexobx = 49, iarr_indexobxprocs = 50, iarr_indexoby = 51, &
& iarr_indexobyprocs = 52, iarr_iobu = 53, iarr_iobuloc = 54, iarr_iobv = 55, &
& iarr_iobvloc = 56, iarr_iobx = 57, iarr_iobxloc = 58, iarr_ioby = 59, &
& iarr_iobyloc = 60, iarr_jobu = 61, iarr_jobuloc = 62, iarr_jobv = 63, &
& iarr_jobvloc = 64, iarr_jobx = 65, iarr_jobxloc = 66, iarr_joby = 67, &
& iarr_jobyloc = 68, iarr_maskatc_int = 69, iarr_mgvars_nc1procs = 70, &
& iarr_mgvars_nc2procs = 71, iarr_mgvars_nr1procs = 72, &
& iarr_mgvars_nr2procs = 73, iarr_ncprocs = 74, iarr_nobuprocs = 75, &
& iarr_nobvprocs = 76, iarr_nobxprocs = 77, iarr_nobyprocs = 78, &
& iarr_nodeatc = 79, iarr_nodeatu = 80, iarr_nodeatuv = 81, &
& iarr_nodeatuw = 82, iarr_nodeatv = 83, iarr_nodeatvw = 84, &
& iarr_node2du = 85, iarr_node2duv = 86, iarr_node2dv = 87, &
& iarr_nosbuprocs = 88, iarr_nosbvprocs = 89, iarr_nrprocs = 90, &
& iarr_nrvbuprocs = 91, iarr_nrvbvprocs = 92, iarr_rlxobcatu = 93, &
& iarr_rlxobcatv = 94, iarr_seapoint = 95, iarr_seapointglb = 96, &
& iarr_soutobv = 97, iarr_soutoby = 98, iarr_westobu = 99, iarr_westobx = 100

!---depths
INTEGER, PARAMETER :: &
& iarr_depmeanatc = 101, iarr_depmeanatu = 102, iarr_depmeanatuv = 103, &
& iarr_depmeanatv = 104, iarr_depmeanglb = 105, iarr_deptotatc = 106, &
& iarr_deptotatc_err = 107, iarr_deptotatc_old = 108, iarr_deptotatu = 109, &
& iarr_deptotatu_old = 110, iarr_deptotatuv = 111, iarr_deptotatv = 112, &
& iarr_deptotatv_old = 113, iarr_dzeta = 114, iarr_zeta = 115, &
& iarr_zeta_old = 116

!---currents
INTEGER, PARAMETER :: &
& iarr_hdvelmag = 117, iarr_hmvelmag = 118, iarr_hmvelmag_hadv_cfl = 119, &
& iarr_hvelmag = 120, iarr_hvelmag_hadv_cfl = 121, iarr_p2dbcgradatu = 122, &
& iarr_p2dbcgradatv = 123, iarr_p3dbcgradatu = 124, iarr_p3dbcgradatv = 125, &
& iarr_udevint = 126, iarr_udfvel = 127, iarr_udvel = 128, &
& iarr_udvel_old = 129, iarr_ufvel = 130, iarr_umpred = 131, iarr_umvel = 132, &
& iarr_umvel_hadv_cfl = 133, iarr_umvel_old = 134, iarr_uvel = 135, &
& iarr_uvel_hadv_cfl = 136, iarr_uvel_old = 137, iarr_vdevint = 138, &
& iarr_vdfvel = 139, iarr_vdvel = 140, iarr_vdvel_old = 141, iarr_vel2d = 142, &
& iarr_vel3d = 143, iarr_vfvel = 144, iarr_vmpred = 145, iarr_vmvel = 146, &
& iarr_vmvel_hadv_cfl = 147, iarr_vmvel_old = 148, iarr_vvel = 149, &
& iarr_vvel_hadv_cfl = 150, iarr_vvel_old = 151, iarr_wfvel = 152, &
& iarr_wphys = 153, iarr_wvel = 154, iarr_wvel_vadv_cfl = 155

!---density
INTEGER, PARAMETER :: &
& iarr_beta_sal = 156, iarr_beta_temp = 157, iarr_dens = 158, iarr_sal = 159, &
& iarr_temp = 160

!---diffusion coefficients
INTEGER, PARAMETER :: &
& iarr_hdifcoef2datc = 161, iarr_hdifcoef2d_mom = 162, &
& iarr_hdifcoef2d_scal = 163, iarr_hdifcoef2datuv = 164, &
& iarr_hdifcoef3datc = 165, iarr_hdifcoef3d_mom = 166, &
& iarr_hdifcoef3d_scal = 167, iarr_hdifcoef3datu = 168, &
& iarr_hdifcoef3datuv = 169, iarr_hdifcoef3datv = 170, &
& iarr_hdifcoef3datw = 171, iarr_kinvisc = 172, iarr_mom_vdif_cfl = 173, &
& iarr_scal_vdif_cfl = 174, iarr_xslopeatu_geo = 175, &
& iarr_xslopeatu_ziso = 176, iarr_xslopeatw_geo = 177, &
& iarr_yslopeatv_geo = 178, iarr_yslopeatv_ziso = 179, &
& iarr_yslopeatw_geo = 180, iarr_vdifcoefmom = 181, iarr_vdifcoefscal = 182, &
& iarr_vdifcoefscal_norot = 183, iarr_vdifcoefscal_rot = 184, &
& iarr_vdifcoeftke = 185, iarr_2D_hdif_cfl = 186

!---turbulence
INTEGER, PARAMETER :: &
& iarr_buofreq2 = 187, iarr_dissip = 188, iarr_ricnum = 189, &
& iarr_shearfreq2 = 190, iarr_tke = 191, iarr_tke_old = 192, iarr_tkezl = 193, &
& iarr_zlmix = 194

!---tides
INTEGER, PARAMETER :: &
& iarr_astro_earth = 195, iarr_fnodal_anal = 196, iarr_fnodal_astro = 197, &
& iarr_fnodal_obc = 198, iarr_fxastro = 199, iarr_fyastro = 200, &
& iarr_index_astro = 201, iarr_index_obc = 202, iarr_ispec_tides = 203, &
& iarr_phase_anal = 204, iarr_phase_astro = 205, iarr_phase_obc = 206, &
& iarr_tidal_spectrum = 207

!---meteo
INTEGER, PARAMETER :: &
& iarr_airtemp = 208, iarr_atmpres = 209, iarr_cloud_cover = 210, &
& iarr_evapminprec = 211, iarr_evaporation = 212, iarr_meteodata = 213, &
& iarr_precipitation = 214, iarr_relhum = 215, iarr_sst = 216, &
& iarr_uwindatc = 217, iarr_uwindatc_old = 218, iarr_vappres_air = 219, &
& iarr_vwindatc = 220, iarr_vwindatc_old = 221, iarr_windatc = 222

!---waves
INTEGER, PARAMETER :: &
& iarr_gxcoordglbwav = 223, iarr_gycoordglbwav = 224, iarr_hbwdissipmag = 225, &
& iarr_hmbwdissipmag = 226, iarr_hmswdissipmag = 227, iarr_hswdissipmag = 228, &
& iarr_maskglbwav = 229, iarr_ubwdissipatc = 230, iarr_umbwdissipatc = 231, &
& iarr_umswdissipatc = 232, iarr_uswdissipatc = 233, iarr_vbwdissipatc = 234, &
& iarr_vmbwdissipatc = 235, iarr_vmswdissipatc = 236, iarr_vswdissipatc = 237, &
& iarr_wavedata = 238, iarr_wavedir = 239, iarr_waveexcurs = 240, &
& iarr_wavefreq = 241, iarr_waveheight = 242, iarr_wavenum = 243, &
& iarr_waveperiod = 244, iarr_wavepres = 245, iarr_wavevel = 246

!---Stokes velocities
INTEGER, PARAMETER :: &
& iarr_hmstokesmag = 247, iarr_hmveltotmag = 248, iarr_hstokesmag = 249, &
& iarr_hveltotmag = 250, iarr_stokessource2du = 251, &
& iarr_stokessource2dv = 252, iarr_udstokesatu = 253, iarr_umstokesatc = 254, &
& iarr_umstokesatu = 255, iarr_umveltot = 256, iarr_ustokesatc = 257, &
& iarr_ustokesatu = 258, iarr_uveltot = 259, iarr_vdstokesatv = 260, &
& iarr_vmstokesatc = 261, iarr_vmstokesatv = 262, iarr_vmveltot = 263, &
& iarr_vstokesatc = 264, iarr_vstokesatv = 265, iarr_vveltot = 266, &
& iarr_wstokesatw = 267

!---optical arrays
INTEGER, PARAMETER :: &
& iarr_optattcoef2 = 268, iarr_qrad = 269, iarr_radiance = 270

!---bottom/surface fluxes
INTEGER, PARAMETER :: &
& iarr_bdragcoefatc = 271, iarr_bfricatu = 272, iarr_bfricatv = 273, &
& iarr_bfricvel = 274, iarr_bfricvel_max = 275, iarr_bfricvel_wav = 276, &
& iarr_bstresatc = 277, iarr_bstresatc_max = 278, iarr_bstresatc_wav = 279, &
& iarr_bstresatu = 280, iarr_bstresatv = 281, iarr_cds = 282, iarr_ces = 283, &
& iarr_chs = 284, iarr_fwave = 285, iarr_qlatent = 286, iarr_qlwave = 287, &
& iarr_qnonsol = 288, iarr_qsensible = 289, iarr_qtot = 290, &
& iarr_sfricatc = 291, iarr_ssalflux = 292, iarr_sstresatc = 293, &
& iarr_ubstresatc = 294, iarr_ubstresatu = 295, iarr_usstresatc = 296, &
& iarr_usstresatu = 297, iarr_vbstresatc = 298, iarr_vbstresatv = 299, &
& iarr_vsstresatc = 300, iarr_vsstresatv = 301, iarr_wavethickatc = 302, &
& iarr_zaroughatc = 303, iarr_zroughatc = 304

!---open boundary forcing
INTEGER, PARAMETER :: &
& iarr_floutobu = 305, iarr_floutobv = 306, iarr_gxslope = 307, &
& iarr_gxslope_amp = 308, iarr_gxslope_pha = 309, iarr_gyslope = 310, &
& iarr_gyslope_amp = 311, iarr_gyslope_pha = 312, iarr_iloczobu = 313, &
& iarr_iloczobv = 314, iarr_indexprof = 315, iarr_index2dobuv = 316, &
& iarr_index2dobxy = 317, iarr_iobc2dtype = 318, iarr_iprofobu = 319, &
& iarr_iprofobv = 320, iarr_iprofobx = 321, iarr_iprofoby = 322, &
& iarr_iqsecobu = 323, iarr_itypobu = 324, iarr_itypobv = 325, &
& iarr_itypobx = 326, iarr_itypoby = 327, iarr_ityp2dobu = 328, &
& iarr_ityp2dobv = 329, iarr_ityp2dobx = 330, iarr_ityp2doby = 331, &
& iarr_jqsecobv = 332, iarr_noprofsd = 333, iarr_no2dobuv = 334, &
& iarr_no2dobxy = 335, iarr_obcdatuv2d = 336, iarr_obcdatxy2d = 337, &
& iarr_obcsalatu = 338, iarr_obcsalatv = 339, iarr_obctmpatu = 340, &
& iarr_obctmpatv = 341, iarr_obc2uvatu = 342, iarr_obc2uvatu_old = 343, &
& iarr_obc2uvatv = 344, iarr_obc2uvatv_old = 345, iarr_obc3uvatu = 346, &
& iarr_obc3uvatv = 347, iarr_prof3dsal = 348, iarr_prof3dtmp = 349, &
& iarr_prof3dveluv = 350, iarr_prof3dvelxy = 351, iarr_return_time = 352, &
& iarr_surdat1d = 353, iarr_udatobu = 354, iarr_udatobu_amp = 355, &
& iarr_udatobu_pha = 356, iarr_udatoby = 357, iarr_udatoby_amp = 358, &
& iarr_udatoby_pha = 359, iarr_vdatobv = 360, iarr_vdatobv_amp = 361, &
& iarr_vdatobv_pha = 362, iarr_vdatobx = 363, iarr_vdatobx_amp = 364, &
& iarr_vdatobx_pha = 365, iarr_zdatobu = 366, iarr_zdatobu_amp = 367, &
& iarr_zdatobu_pha = 368, iarr_zdatobv = 369, iarr_zdatobv_amp = 370, &
& iarr_zdatobv_pha = 371, iarr_zeta_amp = 372, iarr_zeta_pha = 373

!---structures
INTEGER, PARAMETER :: &
& iarr_idry = 374, iarr_indexwbaru = 375, iarr_indexwbaruprocs = 376, &
& iarr_indexwbarv = 377, iarr_indexwbarvprocs = 378, iarr_ithinu = 379, &
& iarr_ithinuloc = 380, iarr_ithinv = 381, iarr_ithinvloc = 382, &
& iarr_iwbaru = 383, iarr_iwbaruloc = 384, iarr_iwbarv = 385, &
& iarr_iwbarvloc = 386, iarr_jdry = 387, iarr_jthinu = 388, &
& iarr_jthinuloc = 389, iarr_jthinv = 390, iarr_jthinvloc = 391, &
& iarr_jwbaru = 392, iarr_jwbaruloc = 393, iarr_jwbarv = 394, &
& iarr_jwbarvloc = 395, iarr_nowbaruprocs = 396, iarr_nowbarvprocs = 397, &
& iarr_oricoefu = 398, iarr_oricoefv = 399, iarr_oriheightu = 400, &
& iarr_oriheightv = 401, iarr_oriwidthu = 402, iarr_oriwidthv = 403, &
& iarr_wbarcoefu = 404, iarr_wbarcoefv = 405, iarr_wbarcrestu = 406, &
& iarr_wbarcrestv = 407, iarr_wbarmodlu = 408, iarr_wbarmodlv = 409, &
& iarr_wbarelossu = 410, iarr_wbarelossv = 411, iarr_mpvcov = 551

!---discharges
INTEGER, PARAMETER :: &
& iarr_disarea = 412, iarr_disdir = 413, iarr_disflag = 414, &
& iarr_dissal = 415, iarr_disspeed = 416, iarr_distmp = 417, &
& iarr_disvol = 418, iarr_idis = 419, iarr_idisloc = 420, &
& iarr_indexdisloc = 421, iarr_indexdisprocs = 422, iarr_jdis = 423, &
& iarr_jdisloc = 424, iarr_kdis = 425, iarr_kdistype = 426, &
& iarr_nodisprocs = 427, iarr_xdiscoord = 428, iarr_ydiscoord = 429, &
& iarr_zdiscoord = 430

!---parameters for parallel mode
INTEGER, PARAMETER :: &
& iarr_idprocs = 431

!---energy equation, enstrophy, vorticity
INTEGER, PARAMETER :: &
& iarr_edens0d = 432, iarr_edens2d = 433, iarr_edens3d = 434, &
& iarr_edissip0d = 435, iarr_edissip2d = 436, iarr_edissip3d = 437, &
& iarr_eflux2du = 438, iarr_eflux2dv = 439, iarr_eflux3du = 440, &
& iarr_eflux3dv = 441, iarr_eflux3dw = 442, iarr_ekin0d = 443, &
& iarr_ekin2d = 444, iarr_ekin3d = 445, iarr_enstr0d = 446, iarr_epot0d = 447, &
& iarr_epot2d = 448, iarr_etot0d = 449, iarr_etot2d = 450, iarr_etot3d = 451, &
& iarr_vortic2d = 452, iarr_vortic3d = 453

!---nesting
INTEGER, PARAMETER :: &
& iarr_inst2dtype = 454, iarr_lbhnstatc = 455, iarr_lbhnstatu = 456, &
& iarr_lbhnstatv = 457, iarr_lbhnstatx = 458, iarr_lbhnstaty = 459, &
& iarr_nohnstatc = 460, iarr_nohnstatu = 461, iarr_nohnstatv = 462, &
& iarr_nohnstatx = 463, iarr_nohnstaty = 464, iarr_nohnstcprocs = 465, &
& iarr_nohnstglbc = 466, iarr_nohnstglbu = 467, iarr_nohnstglbv = 468, &
& iarr_nohnstglbx = 469, iarr_nohnstglby = 470, iarr_nohnstuvprocs = 471, &
& iarr_nohnstxyprocs = 472, iarr_novnst = 473

!---elliptic parameters
INTEGER, PARAMETER :: &
& iarr_ellac2d = 474, iarr_ellac3d = 475, iarr_ellcc2d = 476, &
& iarr_ellcc3d = 477, iarr_ellinc2d = 478, iarr_ellinc3d = 479, &
& iarr_ellip2d = 480, iarr_ellip3d = 481, iarr_ellmaj2d = 482, &
& iarr_ellmaj3d = 483, iarr_ellmin2d = 484, iarr_ellmin3d = 485, &
& iarr_ellpha2d = 486, iarr_ellpha3d = 487

!---relative coordinate arrays
INTEGER, PARAMETER :: &
& iarr_icoordC = 488, iarr_icoordCC = 489, iarr_icoordCU = 490, &
& iarr_icoordCV = 491, iarr_icoordU = 492, iarr_icoordUU = 493, &
& iarr_icoordUY = 494, iarr_icoordV = 495, iarr_icoordVV = 496, &
& iarr_icoordVX = 497, iarr_jcoordC = 498, iarr_jcoordCC = 499, &
& iarr_jcoordCU = 500, iarr_jcoordCV = 501, iarr_jcoordU = 502, &
& iarr_jcoordUU = 503, iarr_jcoordUY = 504, iarr_jcoordV = 505, &
& iarr_jcoordVV = 506, iarr_jcoordVX = 507, iarr_weightsC = 508, &
& iarr_weightsCC = 509, iarr_weightsCU = 510, iarr_weightsCV = 511, &
& iarr_weightsU = 512, iarr_weightsUU = 513, iarr_weightsUY = 514, &
& iarr_weightsV = 515, iarr_weightsVV = 516, iarr_weightsVX = 517

!---(absolute) data coordinate arrays
INTEGER, PARAMETER :: &
& iarr_depmean = 518, iarr_scoord = 519, iarr_scoordatc = 520, &
& iarr_scoordatu = 521, iarr_scoordatv = 522, iarr_scoordatx = 523, &
& iarr_scoordaty = 524, iarr_seamask = 525, iarr_statnames = 526, &
& iarr_subname = 527, iarr_time = 528, iarr_xcoord = 529, iarr_xcoordatc = 530,&
& iarr_xcoordatu = 531, iarr_xcoordatv = 532, iarr_xcoordatx = 533, &
& iarr_xcoordaty = 534, iarr_ycoord = 535, iarr_ycoordatc = 536, &
& iarr_ycoordatu = 537, iarr_ycoordatv = 538, iarr_ycoordatx = 539, &
& iarr_ycoordaty = 540, iarr_zcoord = 541, iarr_zcoordatc = 542, &
& iarr_zcoordatu = 543, iarr_zcoordatv = 544, iarr_zcoordatx = 545, &
& iarr_zcoordaty = 546, iarr_zcoordm = 547, iarr_zetout = 548

!---model parameters
INTEGER, PARAMETER :: &
& iarr_density_ref = 549, iarr_gacc_mean = 550


END MODULE modids
