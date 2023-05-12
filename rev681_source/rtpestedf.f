      subroutine rtpestedf
      
!!     ~ ~ ~ PURPOSE ~ ~ ~
!!     this subroutine computes the daily stream chemical balance
!!     (soluble and sorbed)     

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name          |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    ch_l2(:)      |km            |length of main channel
!!    ch_w(2,:)     |m             |average width of main channel
!!    chpst_conc(:) |mg/(m**3)     |initial chemical concentration in reach
!!    chpst_koc(:)  |m**3/g        |chemical partition coefficient between
!!                                 |water and sediment in reach
!!    chpst_mix(:)  |m2/day        |mixing velocity (diffusion/dispersion) for
!!                                 |chemical in reach
!!    chpst_rea(:)  |1/day         |chemical reaction coefficient in reach
!!    chpst_rsp(:)  |m/day         |resuspension velocity in reach for chemical
!!                                 |sorbed to sediment
!!    chpst_stl(:)  |m/day         |settling velocity in reach for chemical
!!                                 |sorbed to sediment
!!    chpst_vol(:)  |m/day         |chemical volatilization coefficient in 
!!                                 |reach
!!    drift(:)      |kg            |amount of chemical drifting onto main
!!                                 |channel in subbasin
!!    hru_sub(:)    |none          |subbasin number where reach is located
!!    inum1         |none          |reach number
!!    inum2         |none          |inflow hydrograph storage location number
!!    rchdep        |m             |depth of flow on day
!!    rchwtr        |m^3 H2O       |water stored in reach at beginning of day
!!    rnum1         |none          |fraction of overland flow
!!    rtwtr         |m^3 H2O       |water leaving reach on day
!!    sedpst_act(:) |m             |depth of active sediment layer in reach for
!!                                 |chemical
!!    sedpst_bry(:) |m/day         |chemical burial velocity in river bed
!!                                 |sediment
!!    sedpst_conc(:)|mg/(m**3)     |inital chemical concentration in river bed
!!                                 |sediment
!!    sedpst_rea(:) |1/day         |chemical reaction coefficient in river bed
!!                                 |sediment
!!    varoute(11,:) |mg pst        |chemical in solution
!!    varoute(12,:) |mg pst        |chemical sorbed to sediment
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    bury        |mg pst        |loss of chemical from active sediment layer
!!                               |by burial
!!    difus       |mg pst        |diffusion of chemical from sediment to reach
!!    reactb      |mg pst        |amount of chemical in sediment that is lost
!!                               |through reactions
!!    reactw      |mg pst        |amount of chemical in reach that is lost
!!                               |through reactions
!!    resuspst    |mg pst        |amount of chemical moving from sediment to
!!                               |reach due to resuspension
!!    setlpst     |mg pst        |amount of chemical moving from water to
!!                               |sediment due to settling
!!    solpesto    |mg pst/m^3    |soluble chemical concentration in outflow
!!                               |on day
!!    sorpesto    |mg pst/m^3    |sorbed chemical concentration in outflow
!!                               |on day
!!    volatpst    |mg pst        |amount of chemical in reach lost by
!!                               |volatilization
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    bedvol      |m^3           |volume of river bed sediment
!!    chpstmass   |mg pst        |mass of chemical in reach
!!    depth       |m             |depth of water in reach
!!    fd2         |
!!    frsol       |none          |fraction of chemical in reach that is soluble
!!    frsrb       |none          |fraction of chemical in reach that is sorbed
!!    jrch        |none          |reach number
!!    pstin       |mg pst        |total chemical transported into reach
!!                               |during time step
!!    sedcon      |kg/m^3         |sediment concentration
!!    sedpstmass  |mg pst        |mass of chemical in bed sediment
!!    solpstin    |mg pst        |soluble chemical entering reach during 
!!                               |time step
!!    sorpstin    |mg pst        |sorbed chemical entering reach during
!!                               |time step
!!    tday        |days          |flow duration
!!    wtrin       |m^3 H2O       |volume of water entering reach during time
!!                               |step
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    Intrinsic: abs

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

      use parm

      integer :: jrch
      real*8 :: solpstin, sorpstin, pstin, depth, chpstmass, frsol, frsrb
      real*8 :: sedpstmass, bedvol, fd2, wtrin, solmax, sedcon, tday
!! LEP EDF 2021 local variables for dynamic chem-sed routines
      real*8 :: sedinflowst, sedinflow, sedstini, seddep, sedresus, sednet
      real*8 :: bedcumul, bedcumulini, sedoutflow, porosity, density
      real*8 :: seda, sedb, bursedpstmass, kd, chemssconc
      real*8 :: chemresus, chemdep, chembury, chemunbury, sedpst_conca
      real*8 :: bedarea
      real :: apprxzero
      
      apprxzero = 1.e-16

      jrch = 0
      jrch = inum1

!! initialize depth of water for chemical calculations
      depth = 0.
      if (rchdep < 0.1) then
        depth = .1
      else
        depth = rchdep
      endif

!! calculate volume of active river bed sediment layer
      bedarea = 0.
      bedvol = 0.
      bedarea = ch_w(2,jrch) * ch_l2(jrch) * 1000. 
      bedvol = bedarea * sedpst_act(jrch)

!! calculate volume of water entering reach
      wtrin = 0.
      wtrin = varoute(2,inum2) * (1. - rnum1)

!! LEP EDF 2021 prepare sediment fluxes needed for sorbed chemical processes      
      sedinflowst = sedinorg * 1000. !kg incoming (t) plus stored suspended sed from previous timestep (t-1)
      sedinflow = (varoute(3,inum2) * (1. - rnum1)) * 1000. !kg incoming routed sediment from upstream this timestep
      sedstini = 0.
      if (sedinflowst - sedinflow > 0.) then
        sedstini = sedinflowst - sedinflow ! starting (t-1) suspended sediment before inflow, dep, resus
      endif
      seddep = rchdy(57,jrch) * 1000. ! sediment deposited this time step
      sedresus = rchdy(56,jrch) * 1000. ! sediment resuspended this time step
      sednet = seddep - sedresus ! net sediment transfer
      bedcumul = depch(jrch) * 1000. !kg cumulative bed sediment after inflow, dep, resus
      bedcumulini = 0.
      if (bedcumul - seddep + sedresus > 0.) then
        bedcumulini = bedcumul - seddep + sedresus ! starting cumulative bed sediment
      endif
      sedoutflow = sedrch * 1000. ! sediment outflow this time step
      !! Two layer benthic sediment model, Active and Buried
      porosity = 0.5 !! default VVWM porosity fraction
      density = 1350. !! default VVWM benthic sediment density kg/m3
      seda = bedvol* porosity * density !! active layer sediment mass kg
      sedb = 0.
      if (bedcumulini - seda > 0.) then
        sedb = bedcumulini - seda !! buried layer sediment mass kg
      endif
      
!! chemical transported into reach during day
      solpstin = 0.
      sorpstin = 0.
      pstin = 0.
      solpstin = varoute(11,inum2) * (1. - rnum1)
      sorpstin = varoute(12,inum2) * (1. - rnum1)
      pstin = solpstin + sorpstin

!! add chemical drifting from HRUs in subbasin to reach
!      if (rtwtr > 0.) then
!        pstin = pstin + (drift(jrch) * 1.e6)
!      else
!        sedpst_conc(jrch) = sedpst_conc(jrch) + drift(jrch) * 1.e6 /    &
!     &                                                            bedvol
!      endif
 
      !! calculate mass of chemical in reach
      chpstmass = 0.
      chpstmass = pstin + chpst_conc(jrch) * rchwtr
      
      !! calculate mass of chemical in active layer bed sediment
      sedpstmass = 0.
      sedpstmass = sedpst_conc(jrch) * bedvol
      bursedpstmass = sedpst_concb(jrch)*sedb

      if (chpstmass + sedpstmass < 1.e-6) then
        chpst_conc(jrch) = 0.
        sedpst_conc(jrch) = 0.
      end if
      if (chpstmass + sedpstmass < 1.e-6) return

!!in-stream processes
      if (rtwtr / 86400. > 0.002) then
!! LEP EDF 2021 calculate suspended sediment concentration after sediment inflow but before dep/resus/outflow        
        sedcon = 0.
        sedcon = sedinflowst /(rchwtr + wtrin)  !! kg/m3

        !! calculate fraction of soluble and sorbed chemical
        frsol = 0.
        frsrb = 0.
!! LEP EDF 2023 rtpest.f expects a chpst_koc to be a Kd(=Koc * foc) in units of m3/g
!! so we'll use the same and assume the input has been pre-multiplied by default 0.04 foc in soil and in m3/g
        kd = chpst_koc(jrch) * 1000. !! convert m3/g to m3/kg
        
!! LEP EDF 2021 repartition even if there is no incoming chem load, water/sed ratio will change chem equilib
          if (chpst_koc(jrch) > 0.) then
            frsol = 1. / (1. + kd * sedcon)
            fsedsol = 1. / (1. + kd * seda/bedvol) !! active bed aqueous fraction
          else
            frsol = 1. !! default to 100% dissovled if Koc not set
            fsedsol = 1.
          end if
          frsrb = 1. - frsol
        
        chemssconc = chpstmass*frsrb/sedinflowst !! sorbed chemical concentration in sus sed mg/kg ss
        
        !! calculate flow duration
        !! tday = 0.
        !! tday = rttime / 24.0
        !! if (tday > 1.0) tday = 1.0
         tday = 1.0 ! assume daily timestep 

        !! calculate amount of chemical that undergoes chemical or
        !! biological degradation on day in reach
        !! MFW, 3/12/12: modify decay to be 1st order
        !! reactw = chpst_rea(jrch) * chpstmass * tday
        reactw = chpstmass - (chpstmass * EXP(-1. * chpst_rea(jrch)     
     &           * tday))
        if (reactw < apprxzero) then
            reactw = 0.
        end if
        
        !! calculate amount of chemical that volatilizes from reach
        volatpst = chpst_vol(jrch) * frsol * chpstmass * tday / depth
        if (volatpst > frsol * chpstmass) then
          volatpst = frsol * chpstmass    
        end if
        if (volatpst < apprxzero) then
          volatpst = 0.
        end if

! LEP EDF 2021 dynamic sorbed chemical processes coupled to actual sediment flux
        chemresus = 0.
        chemdep = 0.
        chembury = 0.
        chemunbury = 0.
        sedpst_conca = sedpst_conc(jrch)/(porosity * density) ! convert active chem concentration from volume to sed mass basis
        If (sednet > seda) then !! all of the active layer is replaced with incoming chem conc and the previous active + remaining dep is buried
          chemdep = chemssconc*sednet
          chembury = chemssconc*(sednet-seda) + sedpstmass
          chemunbury = 0.
          chemresus = 0.

        else if (sednet > 0.) then !! there was a deposit but not complete turnover of active layer
          chemdep = chemssconc*sednet
          chembury = sedpst_conca*sednet !! displaced active layer chem moves to buried layer
          chemunbury = 0.
          chemresus = 0.

        else !! net resuspension
          sednet = -1. * sednet
          if (sednet > (seda + sedb)) then !! there was net erosion greater then cumulative benthic deposits thus far
            chemdep = 0.
            chembury = 0.
            chemresus = sedpstmass + sedpst_concb(jrch)*sedb ! all chem goes back to channel
            chemunbury = sedpst_concb(jrch)*sedb ! all chem is unburied
          else if (sednet > seda) then !! there was net erosion that was greater than active layer but not all of buried
            chemdep = 0.
            chembury = 0.
            chemunbury = sedpst_concb(jrch)*sednet !! an amount equal to net erosion is unburied, refills the active layer, and remainder added to resuspended
            chemresus = sedpstmass + sedpst_concb(jrch)*(sednet-seda) !!all active plus some buried chem reenters channel
          else if (sednet > 0 ) then !! there was net degradation that only partially resuspended active layer
            chemdep = 0.
            chembury = 0.
            chemunbury = sedpst_concb(jrch)*(sednet) !! an amount equal to net erosion is unburied to fill up active lost to resuspension
            chemresus = sedpst_conca*sednet !! net resus chem at active layer concentration
          else !! no net sediment movement
            chemdep = 0.
            chembury = 0.
            chemunbury = 0.
            chemresus = 0.
          end if
        end if

        !! don't let values have more than 3 digit exponent for fix format output
        if (chemdep < apprxzero) then
            chemdep = 0.
        end if
        if (chemresus < apprxzero) then
            chemresus = 0.
        end if
        if (chembury < apprxzero) then
            chembury = 0.
        end if
        
! LEP EDF 2021 fix this to occur between dissolved states in reach and active sed layer
        !! calculate diffusion of chemical between reach and sediment
        !! difus = D*(aqu. conc. reach - aqu. conc. sed)/dx*dt*bedarea
        difus = chpst_mix(jrch) * ((frsol * chpstmass / (rchwtr + wtrin)) - 
     &   (fsedsol * sedpst_conc(jrch))) * tday / depth * bedarea
        if (abs(difus) < apprxzero) then
            difus = 0.
        end if
        
        
! LEP EDF 2021 Now update the mass balance of the channel chemical one time
        loss = reactw + volatpst - difus - chemresus + chemdep
        if (loss > chpstmass) then
          chpstmass = 0.
        else
          chpstmass = chpstmass - loss
        end if
        
        !!update variables for output.rch
        setlpst = chemdep
        
        resuspst = chemresus
        
        bury = chembury
        
        
! LEP EDF 2021 allow the water concentration to exceed solubility 
        !! verify that water concentration is at or below solubility
        !!solmax = 0.
        !!solmax = pest_sol * (rchwtr + wtrin)
        !!if (solmax < chpstmass * frsol) then
        !! sedpstmass = sedpstmass + (chpstmass * frsol - solmax)
        !! chpstmass = chpstmass - (chpstmass * frsol - solmax)
        !!end if
        
      else   
!!insignificant flow
        sedpstmass = sedpstmass + chpstmass
        chpstmass = 0.
      end if

!! LEP EDF 2021 update mass balances of sediment layers
        
      !! calculate loss of chemical from bed sediments by reaction, active layer only
      tday = 1.
      reactb = sedpstmass - (sedpstmass * EXP(-1. * sedpst_rea(jrch)     
     &           * tday))
      if (reactb < apprxzero) then
          reactb = 0.
      end if
      
      loss = reactb-chemdep+chemresus-chemunbury+chembury+difus
      if (loss > sedpstmass) then
        sedpstmass = 0.
      else
        sedpstmass = sedpstmass - loss
      end if

      !! buried benthic chemical mass update
      loss = chemunbury-chembury
      if (loss > bursedpstmass) then
        bursedpstmass = 0.
      else
        bursedpstmass = bursedpstmass - loss
      end if

!! calculate chemical concentrations at end of day
      chpst_conc(jrch) = 0.
      sedpst_conc(jrch) = 0.
      sedpst_concb(jrch) = 0.
      if (rchwtr + wtrin > 1.e-6) then
        chpst_conc(jrch) = chpstmass / (rchwtr + wtrin)
      else
        sedpstmass = sedpstmass + chpstmass
      end if
      sedpst_conc(jrch) = sedpstmass / bedvol
      sedpst_concb(jrch) = bursedpstmass / (bedcumul - seda)

!! calculate amount of chemical transported out of reach
!! LEP EDF 2021 recalculate dissolved/sorbed partitioning after sediment processes
      sedcon = sedoutflow / rtwtr !! kg / m3
      
      if (chpst_koc(jrch) > 0.) then
          frsol = 1. / (1. + kd * sedcon)
      else
          frsol = 1. !! default to 100% dissovled if Koc not set
      end if
      frsrb = 1. - frsol
        
      if (rtwtr / 86400. > 0.002) then             !Claire, corrected to match line 151
        solpesto = chpst_conc(jrch) * frsol
        sorpesto = chpst_conc(jrch) * frsrb
      else
        solpesto = 0.
        sorpesto = 0.
      end if
      
!! don't let the exponent exceed 2 digits for fixed format output
      if (solpesto < apprxzero) then
          solpesto = 0.
      end if
      if (sorpesto < apprxzero) then
          sorpesto = 0.
      end if
      
      return
      end