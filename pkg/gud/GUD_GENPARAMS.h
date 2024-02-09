C This file contains cog templates.
C Code between template and end marker is autogenerated.
C Add new parameters to params.py
C
CCOG[[[cog import cog; from params import *]]]
CCOG[[[end]]] (checksum: d41d8cd98f00b204e9800998ecf8427e)

#ifdef ALLOW_GUD

CBOP
C     !ROUTINE: GUD_GENPARAMS.h
C     !INTERFACE:
C #include GUD_GENPARAMS.h

C     !DESCRIPTION:
C Contains run-time parameters for the gud package
C
C Requires: GUD_SIZE.h

CCOG[[[cog
CCOGfor name,(coms,conds,conde) in coms.items():
CCOG    cog.out(conds)
CCOG    for sfx,decl in coms.items():
CCOG        cog.out('''
CCOG      COMMON /{name}/
CCOG'''.format(name=name+sfx)[1:])
CCOG        cog.out(',\n'.join('''
CCOG{pre}     &    {param}
CCOG{post}'''.format(param=param, pre=pre, post=post)[1:-1] for tp,dims,param,pre,post in decl))
CCOG        for tp,dims,param,pre,post in decl: cog.out('''
CCOG{pre}      {tp} {param}{dims}
CCOG{post}'''.format(param=param, tp=tp, dims=dims, pre=pre, post=post)[:-1])
CCOG        cog.out('\n')
CCOG    cog.out(conde)
CCOG    cog.out('\n')
CCOG]]]
      COMMON /GUD_CONSTANTS_r/
     &    rad2deg
      _RL rad2deg

#ifdef GUD_ALLOW_CARBON
      COMMON /CARBON_CONSTANTS_r/
     &    Pa2Atm,
     &    ptr2mol,
     &    sca1,
     &    sca2,
     &    sca3,
     &    sca4,
     &    sox1,
     &    sox2,
     &    sox3,
     &    sox4,
     &    oA0,
     &    oA1,
     &    oA2,
     &    oA3,
     &    oA4,
     &    oA5,
     &    oB0,
     &    oB1,
     &    oB2,
     &    oB3,
     &    oC0
      _RL Pa2Atm
      _RL ptr2mol
      _RL sca1
      _RL sca2
      _RL sca3
      _RL sca4
      _RL sox1
      _RL sox2
      _RL sox3
      _RL sox4
      _RL oA0
      _RL oA1
      _RL oA2
      _RL oA3
      _RL oA4
      _RL oA5
      _RL oB0
      _RL oB1
      _RL oB2
      _RL oB3
      _RL oC0
#endif

      COMMON /GUD_PARAMS_l/
     &    gud_linFSConserve,
     &    gud_read_phos
      LOGICAL gud_linFSConserve
      LOGICAL gud_read_phos
      COMMON /GUD_PARAMS_i/
     &    gud_seed,
     &    iDEBUG,
     &    jDEBUG,
     &    kDEBUG
      INTEGER gud_seed
      INTEGER iDEBUG
      INTEGER jDEBUG
      INTEGER kDEBUG
      COMMON /GUD_PARAMS_r/
     &    phymin,
     &    katten_w,
     &    katten_chl,
     &    parfrac,
     &    parconv,
     &    tempnorm,
     &    TempAeArr,
     &    TemprefArr,
     &    TempCoeffArr,
     &    alpfe,
     &    scav,
     &    ligand_tot,
     &    ligand_stab,
     &    freefemax,
     &    scav_rat,
     &    scav_inter,
     &    scav_exp,
     &    scav_R_POPPOC,
     &    depthfesed,
     &    fesedflux,
     &    fesedflux_pcm,
     &    R_CP_fesed,
     &    Knita,
     &    Knitb,
     &    PAR_oxi,
     &    Kdoc,
     &    Kdop,
     &    Kdon,
     &    KdoFe,
     &    KPOC,
     &    KPOP,
     &    KPON,
     &    KPOFe,
     &    KPOSi,
     &    wC_sink,
     &    wP_sink,
     &    wN_sink,
     &    wFe_sink,
     &    wSi_sink,
     &    wPIC_sink,
     &    Kdissc,
#ifdef GUD_ALLOW_CARBON
     &    gud_atmos_pCO2,
     &    R_OP,
     &    R_OC,
     &    m3perkg,
     &    surfSaltMinInit,
     &    surfSaltMaxInit,
     &    surfTempMinInit,
     &    surfTempMaxInit,
     &    surfDICMinInit,
     &    surfDICMaxInit,
     &    surfALKMinInit,
     &    surfALKMaxInit,
     &    surfPO4MinInit,
     &    surfPO4MaxInit,
     &    surfSiMinInit,
     &    surfSiMaxInit,
     &    surfSaltMin,
     &    surfSaltMax,
     &    surfTempMin,
     &    surfTempMax,
     &    surfDICMin,
     &    surfDICMax,
     &    surfALKMin,
     &    surfALKMax,
     &    surfPO4Min,
     &    surfPO4Max,
     &    surfSiMin,
     &    surfSiMax,
#endif
     &    diaz_ini_fac,
     &    O2crit,
     &    denit_NP,
     &    denit_NO3,
     &    NO3crit,
     &    PARmin,
     &    chl2nmax,
     &    synthcost,
     &    expPref,
     &    expPalat,
     &    palat_min,
     &    inhib_graz,
     &    inhib_graz_exp,
     &    hillnum,
     &    hollexp,
     &    phygrazmin,
     &    pmaxPON,
     &    pmaxDON,
     &    pcoefO2,
     &    pmaxDIN,
     &    ksatPOM,
     &    ksatDOM,
     &    ksatDIN,
     &    alpha_hydrol,
     &    yod,
     &    yoe,
     &    ynd,
     &    yne,
     &    fnh4,
     &    ynh4,
     &    yonh4,
     &    fno2,
     &    yno2,
     &    yono2,
     &    depthdenit
      _RL phymin
      _RL katten_w
      _RL katten_chl
      _RL parfrac
      _RL parconv
      _RL tempnorm
      _RL TempAeArr
      _RL TemprefArr
      _RL TempCoeffArr
      _RL alpfe
      _RL scav
      _RL ligand_tot
      _RL ligand_stab
      _RL freefemax
      _RL scav_rat
      _RL scav_inter
      _RL scav_exp
      _RL scav_R_POPPOC
      _RL depthfesed
      _RL fesedflux
      _RL fesedflux_pcm
      _RL R_CP_fesed
      _RL Knita
      _RL Knitb
      _RL PAR_oxi
      _RL Kdoc
      _RL Kdop
      _RL Kdon
      _RL KdoFe
      _RL KPOC
      _RL KPOP
      _RL KPON
      _RL KPOFe
      _RL KPOSi
      _RL wC_sink
      _RL wP_sink
      _RL wN_sink
      _RL wFe_sink
      _RL wSi_sink
      _RL wPIC_sink
      _RL Kdissc
#ifdef GUD_ALLOW_CARBON
      _RL gud_atmos_pCO2
      _RL R_OP
      _RL R_OC
      _RL m3perkg
      _RL surfSaltMinInit
      _RL surfSaltMaxInit
      _RL surfTempMinInit
      _RL surfTempMaxInit
      _RL surfDICMinInit
      _RL surfDICMaxInit
      _RL surfALKMinInit
      _RL surfALKMaxInit
      _RL surfPO4MinInit
      _RL surfPO4MaxInit
      _RL surfSiMinInit
      _RL surfSiMaxInit
      _RL surfSaltMin
      _RL surfSaltMax
      _RL surfTempMin
      _RL surfTempMax
      _RL surfDICMin
      _RL surfDICMax
      _RL surfALKMin
      _RL surfALKMax
      _RL surfPO4Min
      _RL surfPO4Max
      _RL surfSiMin
      _RL surfSiMax
#endif
      _RL diaz_ini_fac
      _RL O2crit
      _RL denit_NP
      _RL denit_NO3
      _RL NO3crit
      _RL PARmin
      _RL chl2nmax
      _RL synthcost
      _RL expPref
      _RL expPalat
      _RL palat_min
      _RL inhib_graz
      _RL inhib_graz_exp
      _RL hillnum
      _RL hollexp
      _RL phygrazmin
      _RL pmaxPON
      _RL pmaxDON
      _RL pcoefO2
      _RL pmaxDIN
      _RL ksatPOM
      _RL ksatDOM
      _RL ksatDIN
      _RL alpha_hydrol
      _RL yod
      _RL yoe
      _RL ynd
      _RL yne
      _RL fnh4
      _RL ynh4
      _RL yonh4
      _RL fno2
      _RL yno2
      _RL yono2
      _RL depthdenit

#ifdef GUD_ALLOW_RADTRANS
      COMMON /GUD_RADTRANS_PARAMS_l/
     &    gud_allomSpectra
      LOGICAL gud_allomSpectra
      COMMON /GUD_RADTRANS_PARAMS_c/
     &    gud_waterAbsorbFile,
     &    gud_phytoAbsorbFile,
     &    gud_particleAbsorbFile
      CHARACTER*80 gud_waterAbsorbFile
      CHARACTER*80 gud_phytoAbsorbFile
      CHARACTER*80 gud_particleAbsorbFile
      COMMON /GUD_RADTRANS_PARAMS_i/
     &    gud_selectSolz,
     &    gud_radtrans_kmax
      INTEGER gud_selectSolz
      INTEGER gud_radtrans_kmax
      COMMON /GUD_RADTRANS_PARAMS_r/
     &    gud_refract_water,
     &    gud_rmud_max,
     &    gud_part_size_P,
     &    gud_waveband_edges,
     &    gud_waveband_centers,
     &    gud_radmodThresh,
     &    gud_rmus,
     &    gud_rmuu,
     &    gud_bbmin,
     &    gud_bbw,
     &    gud_lambda_aCDOM,
     &    gud_Sdom,
     &    gud_aCDOM_fac,
     &    gud_aCarCell,
     &    gud_bCarCell,
     &    gud_absorpSlope,
     &    gud_bbbSlope,
     &    gud_scatSwitchSizeLog,
     &    gud_scatSlopeSmall,
     &    gud_scatSlopeLarge
      _RL gud_refract_water
      _RL gud_rmud_max
      _RL gud_part_size_P
      _RL gud_waveband_edges(nlam+1)
      _RL gud_waveband_centers(nlam)
      _RL gud_radmodThresh
      _RL gud_rmus
      _RL gud_rmuu
      _RL gud_bbmin
      _RL gud_bbw
      _RL gud_lambda_aCDOM
      _RL gud_Sdom
      _RL gud_aCDOM_fac
      _RL gud_aCarCell
      _RL gud_bCarCell
      _RL gud_absorpSlope
      _RL gud_bbbSlope
      _RL gud_scatSwitchSizeLog(nlam)
      _RL gud_scatSlopeSmall(nlam)
      _RL gud_scatSlopeLarge(nlam)
#endif

#ifdef GUD_ALLOW_CDOM
      COMMON /GUD_CDOM_PARAMS_r/
     &    fracCDOM,
     &    CDOMdegrd,
     &    CDOMbleach,
     &    PARCDOM,
     &    R_NP_CDOM,
     &    R_FeP_CDOM,
     &    R_CP_CDOM,
     &    CDOMcoeff
      _RL fracCDOM
      _RL CDOMdegrd
      _RL CDOMbleach
      _RL PARCDOM
      _RL R_NP_CDOM
      _RL R_FeP_CDOM
      _RL R_CP_CDOM
      _RL CDOMcoeff
#endif

      COMMON /GUD_DEPENDENT_PARAMS_i/
     &    kMinFeSed,
     &    kMaxFeSed
      INTEGER kMinFeSed
      INTEGER kMaxFeSed

CCOG[[[end]]] (checksum: 5ea99433fbe29ca2bb762c652b038716)

#endif /* ALLOW_GUD */

