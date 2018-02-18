!=================================================================================================================================
! Copyright (c) 2010-2016  Prof. Claus-Dieter Munz 
! This file is part of FLEXI, a high-order accurate framework for numerically solving PDEs with discontinuous Galerkin methods.
! For more information see https://www.flexi-project.org and https://nrg.iag.uni-stuttgart.de/
!
! FLEXI is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
! FLEXI is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
!
! You should have received a copy of the GNU General Public License along with FLEXI. If not, see <http://www.gnu.org/licenses/>.
!=================================================================================================================================
#include "flexi.h"

!==================================================================================================================================
! Contains the routines that set up communicators and control non-blocking communication
!==================================================================================================================================
MODULE MOD_ParticleInit
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------

! Private Part --------------------------------------------------------------------------------------------------------------------

! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE InitParticles
  MODULE PROCEDURE InitParticles
END INTERFACE

INTERFACE FinalizeParticles
  MODULE PROCEDURE FinalizeParticles
END INTERFACE

PUBLIC::InitParticles,FinalizeParticles
!==================================================================================================================================

PUBLIC::DefineParametersParticles
CONTAINS

!==================================================================================================================================
!> Define parameters 
!==================================================================================================================================
SUBROUTINE DefineParametersParticles()
! MODULES
USE MOD_ReadInTools! ,ONLY: prms,addStrListEntry
USE MOD_Particle_Vars
USE Mod_Particle_Globals

IMPLICIT NONE
!
CHARACTER(32)         :: hilf, hilf2
INTEGER               :: iSpec,iInit,nSpeciesInit
!==================================================================================================================================
CALL prms%CreateIntOption('Part-maxParticleNumber', "Max number of particles in part")
CALL prms%CreateIntOption('Part-nBounds', "")
CALL prms%CreateIntOption('Part-nSpecies', "Number of species in part")
CALL prms%CreateIntOption('Part-IterationForMacroVal', "")
CALL prms%CreateIntOption('Particles-NumberOfRandomVectors', "Number of random vectors")
CALL prms%CreateIntOption('RefMappingGuess', "")
CALL prms%CreateIntOption('BezierElevation', "")
CALL prms%CreateIntOption('BezierSampleN', "")
CALL prms%CreateIntOption('Part-nPeriodicVectors', "")
CALL prms%CreateIntOption('PIC-shapefunction-alpha', "")
CALL prms%CreateIntOption('PIC-NbrOfSFdepoFixes', "")
CALL prms%CreateIntOption('PIC-NbrOfSFdepoLayers', "")

CALL prms%CreateLogicalOption('CountNbOfLostParts', "")
CALL prms%CreateLogicalOption('CartesianPeriodic', "")
CALL prms%CreateLogicalOption('DoWriteStateToHDF5', "")
CALL prms%CreateLogicalOption('DoRefMapping', "")
CALL prms%CreateLogicalOption('MeasureTrackTime', "")
CALL prms%CreateLogicalOption('Part-vMPF', "")
CALL prms%CreateLogicalOption('Part-WriteOutputMesh', "")
CALL prms%CreateLogicalOption('Part-WriteMacroValues',"")
CALL prms%CreateLogicalOption('Part-WriteMacroVolumeValues',"")
CALL prms%CreateLogicalOption('Part-WriteMacroSurfaceValues',"")
CALL prms%CreateLogicalOption('Part-WriteFieldsToVTK',"")
CALL prms%CreateLogicalOption('PIC-SFResampleAnalyzeSurfCollis', "")
CALL prms%CreateLogicalOption('PIC-shapefunction-equi', "")
CALL prms%CreateLogicalOption('printMPINeighborWarnings', "")
CALL prms%CreateLogicalOption('PrintSFDepoWarnings', "")

CALL prms%CreateRealOption('Part-DelayTime', "")
CALL prms%CreateRealOption('Particles-ManualTimeStep', "")
CALL prms%CreateRealOption('Part-SafetyFactor', "")
CALL prms%CreateRealOption('Particles-HaloEpsVelo', "")
CALL prms%CreateRealOption('RefMappingEps', "")
CALL prms%CreateRealOption('BezierEpsilonBilinear', "")
CALL prms%CreateRealOption('PIC-shapefunction-radius', "")

CALL prms%CreateRealArrayOption('Part-FIBGMdeltas', "")
CALL prms%CreateRealArrayOption('Part-FactorFIBGM', "")

! No clue how to set this to read it at runtime
nSpeciesInit = 1!GETINT('Part-nSpecies','1')
ALLOCATE(Species(1:nSpeciesInit))

DO iSpec = 1, 1!nSpeciesInit
  WRITE(UNIT=hilf,FMT='(I2)') iSpec
  
  CALL prms%CreateIntOption('Part-Species'//TRIM(ADJUSTL(hilf))//'-nInits',"")
  
  Species(iSpec)%NumberOfInits         = 0!GETINT('Part-Species'//TRIM(ADJUSTL(hilf))//'-nInits','0')
  ALLOCATE(Species(iSpec)%Init(0:Species(iSpec)%NumberOfInits)) 
  
  DO iInit = 0, Species(iSpec)%NumberOfInits
    ! set help characters
    IF(iInit.EQ.0)THEN
      hilf2=TRIM(ADJUSTL(hilf))
    ELSE ! iInit >0
      WRITE(UNIT=hilf2,FMT='(I2)') iInit
      hilf2=TRIM(ADJUSTL(hilf))//'-Init'//TRIM(ADJUSTL(hilf2))
    END IF ! iInit
  
    ! Define GETINT
    CALL prms%CreateIntOption('Part-Species'//TRIM(ADJUSTL(hilf))//'-nInits',"")

    CALL prms%CreateLogicalOption('Part-Species'//TRIM(ADJUSTL(hilf2))//'-UseForInit',"")
    CALL prms%CreateLogicalOption('Part-Species'//TRIM(ADJUSTL(hilf2))//'-UseForEmission',"")

    CALL prms%CreateStringOption('Part-Species'//TRIM(ADJUSTL(hilf2))//'-SpaceIC',"")
    CALL prms%CreateStringOption('Part-Species'//TRIM(ADJUSTL(hilf2))//'-velocityDistribution',"")

!    CALL prms%CreateRealOption('Part-Species'//TRIM(ADJUSTL(hilf2))//'-InflowRiseTime',"")
    CALL prms%CreateRealOption('Part-Species'//TRIM(ADJUSTL(hilf2))//'-PartDensity',"")
    CALL prms%CreateRealOption('Part-Species'//TRIM(ADJUSTL(hilf2))//'-ParticleEmission',"")
    CALL prms%CreateRealOption('Part-Species'//TRIM(hilf2)//'-VeloIC',"")
    CALL prms%CreateRealOption('Part-Species'//TRIM(hilf2)//'-MassIC',"")

    CALL prms%CreateIntOption('Part-Species'//TRIM(ADJUSTL(hilf2))//'-initialParticleNumber',"")
    CALL prms%CreateIntOption('Part-Species'//TRIM(ADJUSTL(hilf2))//'-ParticleEmissionType',"")
    
    CALL prms%CreateRealArrayOption('Part-Species'//TRIM(hilf2)//'-VeloVecIC',"")

  ENDDO
ENDDO
!CALL prms%CreateRealOption('Part-Species1-initialParticleNumber',"")

END SUBROUTINE DefineParametersParticles

!===================================================================================================================================
! Glue Subroutine for particle initialization 
!===================================================================================================================================
SUBROUTINE InitParticles()
! MODULES
USE MOD_Globals!,                   ONLY: MPIRoot,UNIT_STDOUT
USE Mod_Particle_Globals
USE MOD_ReadInTools
USE MOD_IO_HDF5,                    ONLY: AddToElemData,ElementOut
USE MOD_Mesh_Vars,                  ONLY: nElems
USE MOD_LoadBalance_Vars,           ONLY: nPartsPerElem
USE MOD_Particle_Vars,              ONLY: ParticlesInitIsDone,nSpecies!,WriteMacroSurfaceValues
USE MOD_part_emission,              ONLY: InitializeParticleEmission!, InitializeParticleSurfaceflux
! partns
!USE MOD_DSMC_Analyze,               ONLY: InitHODSMC
!USE MOD_DSMC_Init,                  ONLY: InitDSMC

!USE MOD_LD_Init,                    ONLY: InitLD
!USE MOD_LD_Vars,                    ONLY: useLD

! partns
!USE MOD_DSMC_Vars,                  ONLY: useDSMC, DSMC, DSMC_HOSolution,HODSMC
USE MOD_Mesh_Vars,                  ONLY: nElems
!USE MOD_InitializeBackgroundField,  ONLY: InitializeBackgroundField
!USE MOD_PICInterpolation_Vars,      ONLY: useBGField
!USE MOD_Particle_Boundary_Sampling, ONLY: InitParticleBoundarySampling
#ifdef MPI
USE MOD_Particle_MPI,               ONLY: InitParticleCommSize
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
IF(ParticlesInitIsDone)THEN
   SWRITE(*,*) "InitParticles already called."
   RETURN
END IF
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT PARTICLES ...'

IF(.NOT.ALLOCATED(nPartsPerElem))THEN
  ALLOCATE(nPartsPerElem(1:nElems))
  nPartsPerElem=0
  CALL AddToElemData(ElementOut,'nPartsPerElem',LongIntArray=nPartsPerElem(:))
END IF

Call InitParticleGlobals

CALL InitializeVariables()
SWRITE(UNIT_stdOut,'(A)') ' INIT TEST ...'

CALL InitializeParticleEmission()
!CALL InitializeParticleSurfaceflux()

! Initialize volume sampling
!IF(useDSMC .OR. WriteMacroVolumeValues) THEN
!! definition of DSMC sampling values
!  DSMC%SampNum = 0
!  HODSMC%SampleType = TRIM(GETSTR('DSMC-HOSampling-Type','cell_mean'))
!  IF (TRIM(HODSMC%SampleType).EQ.'cell_mean') THEN
!    HODSMC%nOutputDSMC = 1
!    SWRITE(*,*) 'DSMCHO output order is set to 1 for sampling type cell_mean!'
!    ALLOCATE(DSMC_HOSolution(1:10,1,1,1,1:nElems,1:nSpecies))         
!  ELSE
!    HODSMC%nOutputDSMC = GETINT('Particles-DSMC-OutputOrder','1')
!    ALLOCATE(DSMC_HOSolution(1:11,0:HODSMC%nOutputDSMC,0:HODSMC%nOutputDSMC,0:HODSMC%nOutputDSMC,1:nElems,1:nSpecies))         
!  END IF     
!  DSMC_HOSolution = 0.0
!  CALL InitHODSMC()
!END IF

! Initialize surface sampling
!IF (WriteMacroSurfaceValues.OR.DSMC%CalcSurfaceVal) THEN
!  CALL InitParticleBoundarySampling()
!END IF

!IF (useDSMC) THEN
!  CALL  InitDSMC()
!  IF (useLD) CALL InitLD
!ELSE IF (WriteMacroVolumeValues.OR.WriteMacroSurfaceValues) THEN
!  DSMC%ElectronicModel = .FALSE.
!  DSMC%OutputMeshInit  = .FALSE.
!  DSMC%OutputMeshSamp  = .FALSE.
!END IF

#ifdef MPI
! has to be called AFTER InitializeVariables and InitDSMC 
!CALL InitParticleCommSize()
#endif

ParticlesInitIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT PARTICLES DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitParticles

!===================================================================================================================================
! Initialize the variables first 
!===================================================================================================================================
SUBROUTINE InitializeVariables()
! MODULES
USE MOD_Globals!, ONLY:MPIRoot,UNIT_STDOUT,myRank,nProcessors
!USE MOD_Globals_Vars
USE MOD_ReadInTools
USE MOD_Particle_Vars!, ONLY: 
USE MOD_Particle_Boundary_Vars,ONLY:PartBound,nPartBound!,nAdaptiveBC
!USE MOD_Particle_Mesh_Vars    ,ONLY:NbrOfRegions,RegionBounds
!USE MOD_Mesh_Vars,             ONLY:nElems, BoundaryName,BoundaryType, nBCs
!USE MOD_Particle_Surfaces_Vars,ONLY:BCdata_auxSF
!USE MOD_DSMC_Vars,             ONLY:useDSMC, DSMC
USE MOD_Particle_Output_Vars,  ONLY:WriteFieldsToVTK, OutputMesh
!USE MOD_part_MPFtools,         ONLY:DefinePolyVec, DefineSplitVec
!USE MOD_PICInterpolation,      ONLY:InitializeInterpolation
!USE MOD_PICInit,               ONLY:InitPIC

! partns - do I need this?
USE MOD_Particle_Mesh,         ONLY:InitFIBGM,MapRegionToElem

USE MOD_Particle_Tracking_Vars,ONLY:DoRefMapping
USE MOD_Particle_MPI_Vars,     ONLY:SafetyFactor,halo_eps_velo,PartMPI
!USE MOD_part_pressure,         ONLY:ParticlePressureIni,ParticlePressureCellIni
USE MOD_TimeDisc_Vars,         ONLY:TEnd
#ifdef MPI
!USE MOD_Particle_MPI,          ONLY: InitEmissionComm
#endif /*MPI*/
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: iSpec, iInit!, iPartBound, iSeed, iCC
!INTEGER               :: SeedSize, iPBC, iBC, iSwaps, iRegions, iExclude
INTEGER               :: ALLOCSTAT
CHARACTER(32)         :: hilf, hilf2!, hilf3
!CHARACTER(200)        :: tmpString
!LOGICAL               :: TrueRandom, PartDens_OnlyInit                                                  !
!INTEGER,ALLOCATABLE   :: iSeeds(:)
!REAL                  :: iRan, aVec, bVec   ! random numbers for random vectors
!REAL                  :: lineVector(3), v_drift_line, A_ins
!INTEGER               :: iVec, MaxNbrOfSpeciesSwaps,iIMDSpec
!LOGICAL               :: exitTrue,IsIMDSpecies
INTEGER               :: dummy_int
!===================================================================================================================================
! Read print flags
! printRandomSeeds = GETLOGICAL('printRandomSeeds','.FALSE.')
! Read basic particle parameter
PDM%maxParticleNumber = GETINT('Part-maxParticleNumber','1')
!#if (PP_TimeDiscMethod==1)||(PP_TimeDiscMethod==2)||(PP_TimeDiscMethod==6)||(PP_TimeDiscMethod>=501 && PP_TimeDiscMethod<=506)

IF(DoRefMapping)THEN
  ALLOCATE(PartPosRef(1:3,PDM%MaxParticleNumber), STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL abort(&
  __STAMP__&
  ,' Cannot allocate partposref!')
  PartPosRef=-888.
END IF

! predefine random vectors
NumRanVec = GETINT('Particles-NumberOfRandomVectors','100000')
!IF ((usevMPF).OR.(useDSMC)) THEN
!  ALLOCATE(RandomVec(NumRanVec, 3))
!  RandomVec = 0
!  DO iVec = 1, NumRanVec  ! calculation of NumRanVec different Vectors
!    CALL RANDOM_NUMBER(iRan)
!    bVec              = 1 - 2*iRan
!    aVec              = SQRT(1 - bVec**2)
!    RandomVec(iVec,1) = bVec
!    CALL RANDOM_NUMBER(iRan)
!    bVec              = Pi *2 * iRan
!    RandomVec(iVec,2) = aVec * COS(bVec)
!    RandomVec(iVec,3) = aVec * SIN(bVec)
!  END DO
!END IF

ALLOCATE(PartState(1:PDM%maxParticleNumber,1:6)       , &
         LastPartPos(1:PDM%maxParticleNumber,1:3)     , &
         Pt(1:PDM%maxParticleNumber,1:3)              , &
         PartSpecies(1:PDM%maxParticleNumber)         , &
         PDM%ParticleInside(1:PDM%maxParticleNumber)  , &
         PDM%nextFreePosition(1:PDM%maxParticleNumber), &
         PDM%dtFracPush(1:PDM%maxParticleNumber)      , &
         PDM%IsNewPart(1:PDM%maxParticleNumber), STAT=ALLOCSTAT)
IF (ALLOCSTAT.NE.0) THEN
  CALL abort(&
__STAMP__&
  ,'ERROR in particle_init.f90: Cannot allocate Particle arrays!')
END IF
PDM%ParticleInside(1:PDM%maxParticleNumber)  = .FALSE.
PDM%dtFracPush(1:PDM%maxParticleNumber)      = .FALSE.
PDM%IsNewPart(1:PDM%maxParticleNumber)       = .FALSE.
LastPartPos(1:PDM%maxParticleNumber,1:3)     = 0.
PartState                                    = 0.
Pt                                           = 0.
PartSpecies                                  = 0
PDM%nextFreePosition(1:PDM%maxParticleNumber)= 0

! Initialize of adsorbed particle tracking 
KeepWallParticles = .FALSE.

nSpecies = GETINT('Part-nSpecies','1')

!! init varibale MPF per particle
!IF (usevMPF) THEN
!  enableParticleMerge = GETLOGICAL('Part-vMPFPartMerge','.FALSE.')
!  IF (enableParticleMerge) THEN
!    vMPFMergePolyOrder = GETINT('Part-vMPFMergePolOrder','2')
!    vMPFMergeCellSplitOrder = GETINT('Part-vMPFCellSplitOrder','15')
!    vMPFMergeParticleTarget = GETINT('Part-vMPFMergeParticleTarget','0')
!    IF (vMPFMergeParticleTarget.EQ.0) WRITE(*,*) 'vMPFMergeParticleTarget equals zero: no merging is performed!'
!    vMPFSplitParticleTarget = GETINT('Part-vMPFSplitParticleTarget','0')
!    IF (vMPFSplitParticleTarget.EQ.0) WRITE(*,*) 'vMPFSplitParticleTarget equals zero: no split is performed!'
!    vMPFMergeParticleIter = GETINT('Part-vMPFMergeParticleIter','100')
!    vMPF_velocityDistribution = TRIM(GETSTR('Part-vMPFvelocityDistribution','OVDR'))
!    vMPF_relativistic = GETLOGICAL('Part-vMPFrelativistic','.FALSE.')
!    IF(vMPF_relativistic.AND.(vMPF_velocityDistribution.EQ.'MBDR')) THEN
!      CALL abort(&
!__STAMP__&
!      ,'Relativistic handling of vMPF is not possible using MBDR velocity distribution!')
!    END IF
!    ALLOCATE(vMPF_SpecNumElem(1:nElems,1:nSpecies))
!  END IF
!  ALLOCATE(PartMPF(1:PDM%maxParticleNumber), STAT=ALLOCSTAT)
!  IF (ALLOCSTAT.NE.0) THEN
!    CALL abort(&
!__STAMP__&
!    ,'ERROR in particle_init.f90: Cannot allocate Particle arrays!')
!  END IF
!END IF


!nSpecies = CNTSTR('Part-Species-SpaceIC')

IF (nSpecies.LE.0) THEN
  CALL abort(&
__STAMP__&
  ,'ERROR: nSpecies .LE. 0:', nSpecies)
END IF
!PartPressureCell = .FALSE.
!ALLOCATE(Species(1:nSpecies))

DO iSpec = 1, nSpecies
  WRITE(UNIT=hilf,FMT='(I2)') iSpec
  
  ! Define GETINT
  CALL prms%CreateIntOption('Part-Species'//TRIM(ADJUSTL(hilf))//'-nInits',"")
  
  Species(iSpec)%NumberOfInits         = GETINT('Part-Species'//TRIM(ADJUSTL(hilf))//'-nInits','0')
!  ALLOCATE(Species(iSpec)%Init(0:Species(iSpec)%NumberOfInits)) 
  DO iInit = 0, Species(iSpec)%NumberOfInits
    ! set help characters
    IF(iInit.EQ.0)THEN
      hilf2=TRIM(ADJUSTL(hilf))
    ELSE ! iInit >0
      WRITE(UNIT=hilf2,FMT='(I2)') iInit
      hilf2=TRIM(ADJUSTL(hilf))//'-Init'//TRIM(ADJUSTL(hilf2))
    END IF ! iInit
    ! get species values // only once
!    IF(iInit.EQ.0)THEN
!      !General Species Values
!      Species(iSpec)%ChargeIC              = GETREAL('Part-Species'//TRIM(hilf2)//'-ChargeIC','0.')
      Species(iSpec)%MassIC                = GETREAL('Part-Species'//TRIM(hilf2)//'-MassIC','0.')
!      Species(iSpec)%MacroParticleFactor   = GETREAL('Part-Species'//TRIM(hilf2)//'-MacroParticleFactor','1.')
!#if (PP_TimeDiscMethod==120) || (PP_TimeDiscMethod==121) || (PP_TimeDiscMethod==122)
!      Species(iSpec)%IsImplicit            = GETLOGICAL('Part-Species'//TRIM(hilf2)//'-IsImplicit','.FALSE.')
!#endif
!    END IF ! iInit
    
    ! Define GETINIT
!    CALL prms%CreateLogicalOption('Part-Species'//TRIM(ADJUSTL(hilf2))//'-UseForInit',"")
!    CALL prms%CreateLogicalOption('Part-Species'//TRIM(ADJUSTL(hilf2))//'-UseForEmission',"")
!    
!    CALL prms%CreateStringOption('Part-Species'//TRIM(ADJUSTL(hilf2))//'-SpaceIC',"")
!    CALL prms%CreateStringOption('Part-Species'//TRIM(ADJUSTL(hilf2))//'-velocityDistribution',"")
!    
!    CALL prms%CreateRealOption('Part-Species'//TRIM(ADJUSTL(hilf2))//'-InflowRiseTime',"")
!    CALL prms%CreateRealOption('Part-Species'//TRIM(ADJUSTL(hilf2))//'-PartDensity',"")
!    CALL prms%CreateRealOption('Part-Species'//TRIM(ADJUSTL(hilf2))//'-ParticleEmission',"")
!    
!    CALL prms%CreateIntOption('Part-Species'//TRIM(ADJUSTL(hilf2))//'-initialParticleNumber',"")
!!    CALL prms%CreateIntOption('Part-Species1-initialParticleNumber',"")
!    CALL prms%CreateIntOption('Part-Species'//TRIM(ADJUSTL(hilf2))//'-ParticleEmissionType',"")
    
    ! get emission and init data
    Species(iSpec)%Init(iInit)%UseForInit            = GETLOGICAL('Part-Species'//TRIM(ADJUSTL(hilf2))//'-UseForInit','.TRUE.')
    Species(iSpec)%Init(iInit)%UseForEmission        = GETLOGICAL('Part-Species'//TRIM(ADJUSTL(hilf2))//'-UseForEmission','.TRUE.')
    Species(iSpec)%Init(iInit)%SpaceIC               = TRIM(GETSTR('Part-Species'//TRIM(ADJUSTL(hilf2))//'-SpaceIC','cuboid'))
    Species(iSpec)%Init(iInit)%velocityDistribution  = TRIM(GETSTR('Part-Species'//TRIM(ADJUSTL(hilf2))//'-velocityDistribution','constant'))
!    IF(TRIM(Species(iSpec)%Init(iInit)%velocityDistribution).EQ.'tangential_constant')THEN
!      Species(iSpec)%Init(iInit)%Rotation        = GETINT('Part-Species'//TRIM(hilf2)//'-rotation','1')
!      Species(iSpec)%Init(iInit)%VelocitySpread  = GETREAL('Part-Species'//TRIM(hilf2)//'-velocityspread','0.')
!      IF(Species(iSpec)%Init(iInit)%VelocitySpread.LT.0. .OR. Species(iSpec)%Init(iInit)%VelocitySpread.GT.1.) CALL abort(&
!__STAMP__&
!          ,' Wrong input parameter for VelocitySpread in [0;1].')
!      Species(iSpec)%Init(iInit)%VelocitySpreadMethod  = GETINT('Part-Species'//TRIM(hilf2)//'-velocityspreadmethod','0')
!    END IF
!    Species(iSpec)%Init(iInit)%InflowRiseTime        = GETREAL('Part-Species'//TRIM(ADJUSTL(hilf2))//'-InflowRiseTime','0.')
    Species(iSpec)%Init(iInit)%initialParticleNumber = GETINT('Part-Species'//TRIM(ADJUSTL(hilf2))//'-initialParticleNumber','0')
!    Species(iSpec)%Init(iInit)%RadiusIC              = GETREAL('Part-Species'//TRIM(hilf2)//'-RadiusIC','1.')
!    Species(iSpec)%Init(iInit)%Radius2IC             = GETREAL('Part-Species'//TRIM(hilf2)//'-Radius2IC','0.')
!    Species(iSpec)%Init(iInit)%RadiusICGyro          = GETREAL('Part-Species'//TRIM(hilf2)//'-RadiusICGyro','1.')
!    Species(iSpec)%Init(iInit)%NormalIC              = GETREALARRAY('Part-Species'//TRIM(hilf2)//'-NormalIC',3,'0. , 0. , 1.')
!    Species(iSpec)%Init(iInit)%BasePointIC           = GETREALARRAY('Part-Species'//TRIM(hilf2)//'-BasePointIC',3,'0. , 0. , 0.')
!    Species(iSpec)%Init(iInit)%BaseVector1IC         = GETREALARRAY('Part-Species'//TRIM(hilf2)//'-BaseVector1IC',3,'1. , 0. , 0.')
!    Species(iSpec)%Init(iInit)%BaseVector2IC         = GETREALARRAY('Part-Species'//TRIM(hilf2)//'-BaseVector2IC',3,'0. , 1. , 0.')
!    Species(iSpec)%Init(iInit)%CuboidHeightIC        = GETREAL('Part-Species'//TRIM(hilf2)//'-CuboidHeightIC','1.')
!    Species(iSpec)%Init(iInit)%CylinderHeightIC      = GETREAL('Part-Species'//TRIM(hilf2)//'-CylinderHeightIC','1.')
!    Species(iSpec)%Init(iInit)%CalcHeightFromDt      = GETLOGICAL('Part-Species'//TRIM(hilf2)//'-CalcHeightFromDt','.FALSE.')
    Species(iSpec)%Init(iInit)%VeloIC                = GETREAL('Part-Species'//TRIM(hilf2)//'-VeloIC','0.')
    Species(iSpec)%Init(iInit)%VeloVecIC             = GETREALARRAY('Part-Species'//TRIM(hilf2)//'-VeloVecIC',3,'0. , 0. , 0.')
!    Species(iSpec)%Init(iInit)%Amplitude             = GETREAL('Part-Species'//TRIM(hilf2)//'-Amplitude','0.01')
!    Species(iSpec)%Init(iInit)%WaveNumber            = GETREAL('Part-Species'//TRIM(hilf2)//'-WaveNumber','2.')
!    Species(iSpec)%Init(iInit)%maxParticleNumberX    = GETINT('Part-Species'//TRIM(hilf2)//'-maxParticleNumber-x','0')
!    Species(iSpec)%Init(iInit)%maxParticleNumberY    = GETINT('Part-Species'//TRIM(hilf2)//'-maxParticleNumber-y','0')
!    Species(iSpec)%Init(iInit)%maxParticleNumberZ    = GETINT('Part-Species'//TRIM(hilf2)//'-maxParticleNumber-z','0')
!    Species(iSpec)%Init(iInit)%Alpha                 = GETREAL('Part-Species'//TRIM(hilf2)//'-Alpha','0.')
!    Species(iSpec)%Init(iInit)%MWTemperatureIC       = GETREAL('Part-Species'//TRIM(hilf2)//'-MWTemperatureIC','0.')
!    Species(iSpec)%Init(iInit)%ConstantPressure      = GETREAL('Part-Species'//TRIM(hilf2)//'-ConstantPressure','0.')
!    Species(iSpec)%Init(iInit)%ConstPressureRelaxFac = GETREAL('Part-Species'//TRIM(hilf2)//'-ConstPressureRelaxFac','1.')
    Species(iSpec)%Init(iInit)%PartDensity           = GETREAL('Part-Species'//TRIM(ADJUSTL(hilf2))//'-PartDensity','0.')
    Species(iSpec)%Init(iInit)%ParticleEmissionType  = GETINT('Part-Species'//TRIM(ADJUSTL(hilf2))//'-ParticleEmissionType','2')
    Species(iSpec)%Init(iInit)%ParticleEmission      = GETREAL('Part-Species'//TRIM(ADJUSTL(hilf2))//'-ParticleEmission','0.')
!    Species(iSpec)%Init(iInit)%NSigma                = GETREAL('Part-Species'//TRIM(hilf2)//'-NSigma','10.')
!    Species(iSpec)%Init(iInit)%NumberOfExcludeRegions= GETINT('Part-Species'//TRIM(hilf2)//'-NumberOfExcludeRegions','0')
!    Species(iSpec)%Init(iInit)%InsertedParticle      = 0
!    Species(iSpec)%Init(iInit)%InsertedParticleSurplus = 0
!    IF(TRIM(Species(iSpec)%Init(iInit)%velocityDistribution).EQ.'maxwell-juettner') THEN
!      Species(iSpec)%Init(iInit)%MJxRatio       = GETREAL('Part-Species'//TRIM(hilf2)//'-MJxRatio','0')
!      Species(iSpec)%Init(iInit)%MJyRatio       = GETREAL('Part-Species'//TRIM(hilf2)//'-MJyRatio','0')
!      Species(iSpec)%Init(iInit)%MJzRatio       = GETREAL('Part-Species'//TRIM(hilf2)//'-MJzRatio','0')
!    END IF
!    IF(TRIM(Species(iSpec)%Init(iInit)%velocityDistribution).EQ.'weibel') THEN
!      Species(iSpec)%Init(iInit)%WeibelVeloPar       = GETREAL('Part-Species'//TRIM(hilf2)//'-WeibelVeloPar','0')
!      Species(iSpec)%Init(iInit)%WeibelVeloPer       = GETREAL('Part-Species'//TRIM(hilf2)//'-WeibelVeloPer','0')
!    END IF
!    IF(TRIM(Species(iSpec)%Init(iInit)%velocityDistribution).EQ. 'OneD-twostreaminstabilty') THEN
!      Species(iSpec)%Init(iInit)%OneDTwoStreamVelo   = GETREAL('Part-Species'//TRIM(hilf2)//'-OneDTwoStreamVelo','0')
!      Species(iSpec)%Init(iInit)%OneDTwoStreamTransRatio = GETREAL('Part-Species'//TRIM(hilf2)//'-OneDTwoStreamTransRatio','0')
!    END IF
!
    !----------- various checks/calculations after read-in of Species(i)%Init(iInit)%-data ----------------------------------!
    !--- Check if Initial ParticleInserting is really used
    IF ( ((Species(iSpec)%Init(iInit)%ParticleEmissionType.EQ.1).OR.(Species(iSpec)%Init(iInit)%ParticleEmissionType.EQ.2)) &
      .AND.(Species(iSpec)%Init(iInit)%UseForInit) ) THEN
      IF ( (Species(iSpec)%Init(iInit)%initialParticleNumber.EQ.0) &
      .AND. (ABS(Species(iSpec)%Init(iInit)%PartDensity).LE.0.) ) THEN
        Species(iSpec)%Init(iInit)%UseForInit=.FALSE.
        SWRITE(*,*) "WARNING: Initial ParticleInserting disabled as neither ParticleNumber"
        SWRITE(*,*) "nor PartDensity detected for Species, Init ", iSpec, iInit
      END IF
    END IF
!    !--- cuboid-/cylinder-height calculation from v and dt
!    IF (.NOT.Species(iSpec)%Init(iInit)%CalcHeightFromDt) THEN
!      IF (TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'cuboid') THEN
!        IF (ALMOSTEQUAL(Species(iSpec)%Init(iInit)%CuboidHeightIC,-1.)) THEN ! flag is initialized with -1, compatibility issue 
!          Species(iSpec)%Init(iInit)%CalcHeightFromDt=.TRUE.                 
!          SWRITE(*,*) "WARNING: Cuboid height will be calculated from v and dt!"
!        END IF
!      ELSE IF (TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'cylinder') THEN
!        IF (ALMOSTEQUAL(Species(iSpec)%Init(iInit)%CylinderHeightIC,-1.)) THEN !flag is initialized with -1, compatibility issue 
!          Species(iSpec)%Init(iInit)%CalcHeightFromDt=.TRUE.                   
!          SWRITE(*,*) "WARNING: Cylinder height will be calculated from v and dt!"
!        END IF
!      END IF
!    END IF
!    IF (Species(iSpec)%Init(iInit)%CalcHeightFromDt) THEN
!      IF ( (Species(iSpec)%Init(iInit)%ParticleEmissionType.NE.1) .AND. (Species(iSpec)%Init(iInit)%ParticleEmissionType.NE.2) ) &
!        CALL abort(&
!__STAMP__&
!          ,' Calculating height from v and dt is only supported for EmiType1 or EmiType2(=default)!')
!      IF ((TRIM(Species(iSpec)%Init(iInit)%SpaceIC).NE.'cuboid') &
!          .AND.(TRIM(Species(iSpec)%Init(iInit)%SpaceIC).NE.'cylinder')) &
!        CALL abort(&
!__STAMP__&
!          ,' Calculating height from v and dt is only supported for cuboid or cylinder!')
!      IF (Species(iSpec)%Init(iInit)%UseForInit) &
!        CALL abort(&
!__STAMP__&
!          ,' Calculating height from v and dt is not supported for initial ParticleInserting!')
!    END IF
!    !--- virtual pre-insertion (vpi) checks and calculations
!    IF ((TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'cuboid_vpi') &
!      .OR.(TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'cylinder_vpi')) THEN
!      IF ( (Species(iSpec)%Init(iInit)%ParticleEmissionType.NE.1) .AND. (Species(iSpec)%Init(iInit)%ParticleEmissionType.NE.2) ) &
!        CALL abort(&
!__STAMP__&
!        ,' Wrong emission-type for virtual Pre-Inserting region!')
!      IF (TRIM(Species(iSpec)%Init(iInit)%velocityDistribution).NE.'maxwell_lpn') &
!        CALL abort(&
!__STAMP__&
!        ,' Only maxwell_lpn is implemened as velocity-distribution for virtual Pre-Inserting region!')
!      IF (Species(iSpec)%Init(iInit)%UseForInit) &
!        CALL abort(&
!__STAMP__&
!          ,' virtual Pre-Inserting is not supported for initial ParticleInserting. Use additional Init!')
!      !-- Virtual Pre-Inserting is used correctly !
!      Species(iSpec)%Init(iInit)%VirtPreInsert = .TRUE.
!      SWRITE(*,*) "Virtual Pre-Inserting is used for Species, Init ", iSpec, iInit
!      IF (Species(iSpec)%Init(iInit)%PartDensity .EQ. 0.) THEN
!        SWRITE(*,*) "WARNING: If VPI-BC is open, a backflow might not be compensated"
!        SWRITE(*,*) "         (use PartDensity instead of ParticleEmission)!"
!      END IF
!      Species(iSpec)%Init(iInit)%vpiDomainType = TRIM(GETSTR('Part-Species'//TRIM(hilf2)//'-vpiDomainType','perpendicular_extrusion'))
!      SELECT CASE ( TRIM(Species(iSpec)%Init(iInit)%vpiDomainType) )
!      CASE ( 'freestream' )
!        IF ( TRIM(Species(iSpec)%Init(iInit)%SpaceIC) .NE. 'cuboid_vpi' ) THEN
!          CALL abort(&
!__STAMP__&
!            ,' Only cuboid_vpi is supported for a freestream vpiDomainType! (Use default vpiDomainType for cylinder.)')
!        ELSE
!          Species(iSpec)%Init(iInit)%vpiBVBuffer(1) = GETLOGICAL('Part-Species'//TRIM(hilf2)//'-vpiBV1BufferNeg','.TRUE.')
!          Species(iSpec)%Init(iInit)%vpiBVBuffer(2) = GETLOGICAL('Part-Species'//TRIM(hilf2)//'-vpiBV1BufferPos','.TRUE.')
!          Species(iSpec)%Init(iInit)%vpiBVBuffer(3) = GETLOGICAL('Part-Species'//TRIM(hilf2)//'-vpiBV2BufferNeg','.TRUE.')
!          Species(iSpec)%Init(iInit)%vpiBVBuffer(4) = GETLOGICAL('Part-Species'//TRIM(hilf2)//'-vpiBV2BufferPos','.TRUE.')
!        END IF
!      CASE ( 'orifice' )
!        Species(iSpec)%Init(iInit)%vpiBVBuffer = .TRUE.
!        IF ( ABS(Species(iSpec)%Init(iInit)%Radius2IC) .GT. 0. ) THEN
!          CALL abort(&
!__STAMP__&
!            ,' Annular orifice is not implemented yet!')
!        END IF
!      CASE ( 'perpendicular_extrusion' )
!        Species(iSpec)%Init(iInit)%vpiBVBuffer = .TRUE. !dummy
!      CASE DEFAULT
!        CALL abort(&
!__STAMP__&
!,'vpiDomainType is not implemented!')
!      END SELECT
!      !--
!    ELSE
!      Species(iSpec)%Init(iInit)%VirtPreInsert = .FALSE.
!    END IF
!    !--- integer check for ParticleEmissionType 2
!    IF((Species(iSpec)%Init(iInit)%ParticleEmissionType.EQ.2).AND. &
!         ((Species(iSpec)%Init(iInit)%ParticleEmission-INT(Species(iSpec)%Init(iInit)%ParticleEmission)).NE.0)) THEN
!       CALL abort(&
!__STAMP__&
!       ,' If ParticleEmissionType = 2 (parts per iteration), ParticleEmission has to be an integer number')
!    END IF
!    !--- flag for cell-based constant-pressure-EmiTypes
!    IF ((Species(iSpec)%Init(iInit)%ParticleEmissionType .EQ. 4).OR. &
!        (Species(iSpec)%Init(iInit)%ParticleEmissionType .EQ. 6)) PartPressureCell = .TRUE.
!    !--- normalize VeloVecIC and NormalIC (and BaseVector 1 & 2 IC for cylinder) for Inits
!    IF (.NOT. ALL(Species(iSpec)%Init(iInit)%VeloVecIC(:).eq.0.)) THEN
!      Species(iSpec)%Init(iInit)%VeloVecIC = Species(iSpec)%Init(iInit)%VeloVecIC            / &
!        SQRT(Species(iSpec)%Init(iInit)%VeloVecIC(1)*Species(iSpec)%Init(iInit)%VeloVecIC(1) + &
!        Species(iSpec)%Init(iInit)%VeloVecIC(2)*Species(iSpec)%Init(iInit)%VeloVecIC(2)      + &
!        Species(iSpec)%Init(iInit)%VeloVecIC(3)*Species(iSpec)%Init(iInit)%VeloVecIC(3))
!    END IF
!    Species(iSpec)%Init(iInit)%NormalIC = Species(iSpec)%Init(iInit)%NormalIC /                 &
!      SQRT(Species(iSpec)%Init(iInit)%NormalIC(1)*Species(iSpec)%Init(iInit)%NormalIC(1) + &
!      Species(iSpec)%Init(iInit)%NormalIC(2)*Species(iSpec)%Init(iInit)%NormalIC(2) + &
!      Species(iSpec)%Init(iInit)%NormalIC(3)*Species(iSpec)%Init(iInit)%NormalIC(3))
!    IF ((TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'cylinder')&
!        .OR.(TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'cylinder_vpi')) THEN
!        Species(iSpec)%Init(iInit)%BaseVector1IC =&
!                  Species(iSpec)%Init(iInit)%RadiusIC * Species(iSpec)%Init(iInit)%BaseVector1IC /     &
!        SQRT(Species(iSpec)%Init(iInit)%BaseVector1IC(1)*Species(iSpec)%Init(iInit)%BaseVector1IC(1) + &
!        Species(iSpec)%Init(iInit)%BaseVector1IC(2)*Species(iSpec)%Init(iInit)%BaseVector1IC(2) + &
!        Species(iSpec)%Init(iInit)%BaseVector1IC(3)*Species(iSpec)%Init(iInit)%BaseVector1IC(3))
!        Species(iSpec)%Init(iInit)%BaseVector2IC =&
!                   Species(iSpec)%Init(iInit)%RadiusIC * Species(iSpec)%Init(iInit)%BaseVector2IC /    &
!        SQRT(Species(iSpec)%Init(iInit)%BaseVector2IC(1)*Species(iSpec)%Init(iInit)%BaseVector2IC(1) + &
!        Species(iSpec)%Init(iInit)%BaseVector2IC(2)*Species(iSpec)%Init(iInit)%BaseVector2IC(2)      + &
!        Species(iSpec)%Init(iInit)%BaseVector2IC(3)*Species(iSpec)%Init(iInit)%BaseVector2IC(3))
!    END IF
!    !--- read stuff for ExcludeRegions and normalize/calculate corresponding vectors
!    IF (Species(iSpec)%Init(iInit)%NumberOfExcludeRegions.GT.0) THEN
!      ALLOCATE(Species(iSpec)%Init(iInit)%ExcludeRegion(1:Species(iSpec)%Init(iInit)%NumberOfExcludeRegions)) 
!      IF (((TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'cuboid') &
!       .OR.(TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'cylinder')) &
!      .OR.((TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'cuboid_vpi') &
!       .OR.(TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'cylinder_vpi'))) THEN
!        DO iExclude=1,Species(iSpec)%Init(iInit)%NumberOfExcludeRegions
!          WRITE(UNIT=hilf3,FMT='(I2)') iExclude
!          hilf3=TRIM(hilf2)//'-ExcludeRegion'//TRIM(hilf3)
!          Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%SpaceIC             &
!            = TRIM(GETSTR('Part-Species'//TRIM(hilf3)//'-SpaceIC','cuboid'))
!          Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%RadiusIC             &
!            = GETREAL('Part-Species'//TRIM(hilf3)//'-RadiusIC','1.')
!          Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%Radius2IC            &
!            = GETREAL('Part-Species'//TRIM(hilf3)//'-Radius2IC','0.')
!          Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%NormalIC             &
!            = GETREALARRAY('Part-Species'//TRIM(hilf3)//'-NormalIC',3,'0. , 0. , 1.')
!          Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BasePointIC          &
!            = GETREALARRAY('Part-Species'//TRIM(hilf3)//'-BasePointIC',3,'0. , 0. , 0.')
!          Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector1IC        &
!            = GETREALARRAY('Part-Species'//TRIM(hilf3)//'-BaseVector1IC',3,'1. , 0. , 0.')
!          Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector2IC        &
!            = GETREALARRAY('Part-Species'//TRIM(hilf3)//'-BaseVector2IC',3,'0. , 1. , 0.')
!          Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%CuboidHeightIC       &
!            = GETREAL('Part-Species'//TRIM(hilf3)//'-CuboidHeightIC','1.')
!          Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%CylinderHeightIC     &
!            = GETREAL('Part-Species'//TRIM(hilf3)//'-CylinderHeightIC','1.')
!          !--normalize and stuff
!          IF ((TRIM(Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%SpaceIC).EQ.'cuboid') .OR. &
!               ((((.NOT.ALMOSTEQUAL(Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector1IC(1),1.) &
!              .OR. .NOT.ALMOSTEQUAL(Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector1IC(2),0.)) &
!              .OR. .NOT.ALMOSTEQUAL(Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector1IC(3),0.)) &
!            .OR. ((.NOT.ALMOSTEQUAL(Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector2IC(1),0.) &
!              .OR. .NOT.ALMOSTEQUAL(Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector2IC(2),1.)) &
!              .OR. .NOT.ALMOSTEQUAL(Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector2IC(3),0.))) &
!            .AND. (((ALMOSTEQUAL(Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%NormalIC(1),0.)) &
!              .AND. (ALMOSTEQUAL(Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%NormalIC(2),0.))) &
!              .AND. (ALMOSTEQUAL(Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%NormalIC(3),1.))))) THEN
!            !-- cuboid; or BV are non-default and NormalIC is default: calc. NormalIC for ExcludeRegions from BV1/2
!            !   (for def. BV and non-def. NormalIC; or all def. or non-def.: Use User-defined NormalIC when ExclRegion is cylinder)
!            Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%NormalIC(1) &
!              = Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector1IC(2) &
!              * Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector2IC(3) &
!              - Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector1IC(3) &
!              * Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector2IC(2)
!            Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%NormalIC(2) &
!              = Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector1IC(3) &
!              * Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector2IC(1) &
!              - Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector1IC(1) &
!              * Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector2IC(3)
!            Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%NormalIC(3) &
!              = Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector1IC(1) &
!              * Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector2IC(2) &
!              - Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector1IC(2) &
!              * Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector2IC(1)
!          ELSE IF ( (TRIM(Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%SpaceIC).NE.'cuboid') .AND. &
!                    (TRIM(Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%SpaceIC).NE.'cylinder') )THEN
!            CALL abort(&
!__STAMP__&
!,'Error in ParticleInit, ExcludeRegions must be cuboid or cylinder!')
!          END IF
!          IF (Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%NormalIC(1)**2 + &
!              Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%NormalIC(2)**2 + &
!              Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%NormalIC(3)**2 .GT. 0.) THEN
!            Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%NormalIC &
!              = Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%NormalIC &
!              / SQRT(Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%NormalIC(1)**2 &
!              + Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%NormalIC(2)**2 &
!              + Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%NormalIC(3)**2)
!            Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%ExcludeBV_lenghts(1) &
!              = SQRT(Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector1IC(1)**2 &
!              + Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector1IC(2)**2 &
!              + Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector1IC(3)**2)
!            Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%ExcludeBV_lenghts(2) &
!              = SQRT(Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector2IC(1)**2 &
!              + Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector2IC(2)**2 &
!              + Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector2IC(3)**2)
!          ELSE
!            CALL abort(&
!__STAMP__&
!,'Error in ParticleInit, NormalIC Vector must not be zero!')
!          END IF
!        END DO !iExclude
!      ELSE
!        CALL abort(&
!__STAMP__&
!,'Error in ParticleInit, ExcludeRegions are currently only implemented for the SpaceIC cuboid(_vpi) or cylinder(_vpi)!')
!      END IF
!    END IF
!    !--- stuff for calculating ParticleEmission/InitialParticleNumber from PartDensity when this value is not used for LD-stuff
!    !                                                                                  (additional not-LD checks might be necassary)
!    PartDens_OnlyInit=.FALSE.
!    IF ((Species(iSpec)%Init(iInit)%PartDensity.GT.0.).AND.(TRIM(Species(iSpec)%Init(iInit)%SpaceIC).NE.'LD_insert')) THEN
!      IF (Species(iSpec)%Init(iInit)%ParticleEmissionType.NE.1) THEN
!        IF ( (Species(iSpec)%Init(iInit)%ParticleEmissionType.EQ.2) .AND. (Species(iSpec)%Init(iInit)%UseForInit) ) THEN
!          PartDens_OnlyInit=.TRUE.
!        ELSE
!          CALL abort(&
!__STAMP__&
!            , 'PartDensity without LD is only supported for EmiType1 or initial ParticleInserting with EmiType1/2!')
!        END IF
!      END IF
!      IF ((TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'cuboid').OR.(TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'cylinder')) THEN
!        IF  ((((TRIM(Species(iSpec)%Init(iInit)%velocityDistribution).EQ.'constant') &
!          .OR.(TRIM(Species(iSpec)%Init(iInit)%velocityDistribution).EQ.'maxwell') ) &
!          .OR.(TRIM(Species(iSpec)%Init(iInit)%velocityDistribution).EQ.'maxwell_lpn') ) &
!          .OR.(TRIM(Species(iSpec)%Init(iInit)%velocityDistribution).EQ.'emmert') ) THEN
!          IF (Species(iSpec)%Init(iInit)%ParticleEmission .GT. 0.) THEN
!            CALL abort(&
!__STAMP__&
!            ,'Either ParticleEmission or PartDensity can be defined for selected emission parameters, not both!')
!          END IF
!          !---calculation of Base-Area and corresponding component of VeloVecIC
!          lineVector(1) = Species(iSpec)%Init(iInit)%BaseVector1IC(2) * Species(iSpec)%Init(iInit)%BaseVector2IC(3) - &
!            Species(iSpec)%Init(iInit)%BaseVector1IC(3) * Species(iSpec)%Init(iInit)%BaseVector2IC(2)
!          lineVector(2) = Species(iSpec)%Init(iInit)%BaseVector1IC(3) * Species(iSpec)%Init(iInit)%BaseVector2IC(1) - &
!            Species(iSpec)%Init(iInit)%BaseVector1IC(1) * Species(iSpec)%Init(iInit)%BaseVector2IC(3)
!          lineVector(3) = Species(iSpec)%Init(iInit)%BaseVector1IC(1) * Species(iSpec)%Init(iInit)%BaseVector2IC(2) - &
!            Species(iSpec)%Init(iInit)%BaseVector1IC(2) * Species(iSpec)%Init(iInit)%BaseVector2IC(1)
!          A_ins = lineVector(1)*lineVector(1) + lineVector(2)*lineVector(2) + lineVector(3)*lineVector(3)
!          IF (A_ins .GT. 0.) THEN
!            A_ins = SQRT(A_ins)
!            lineVector = lineVector / A_ins
!            IF (Species(iSpec)%Init(iInit)%CalcHeightFromDt) THEN
!              v_drift_line = Species(iSpec)%Init(iInit)%VeloIC * &
!                ( Species(iSpec)%Init(iInit)%VeloVecIC(1)*lineVector(1) + Species(iSpec)%Init(iInit)%VeloVecIC(2)*lineVector(2) &
!                + Species(iSpec)%Init(iInit)%VeloVecIC(3)*lineVector(3) ) !lineVector component of drift-velocity
!            ELSE
!              v_drift_line = 0.
!              IF (Species(iSpec)%Init(iInit)%UseForInit) THEN
!                PartDens_OnlyInit=.TRUE.
!              ELSE
!                CALL abort(&
!__STAMP__&
!                  ,'PartDensity without LD is only supported for CalcHeightFromDt, vpi, or initial ParticleInserting!')
!              END IF
!            END IF
!            IF ( TRIM(Species(iSpec)%Init(iInit)%SpaceIC) .EQ. 'cylinder' ) THEN
!              A_ins = Pi * (Species(iSpec)%Init(iInit)%RadiusIC**2-Species(iSpec)%Init(iInit)%Radius2IC**2)
!            END IF
!            !---calculation of particle flow (macroparticles/s) through boundary
!            IF (.NOT.PartDens_OnlyInit) THEN
!              Species(iSpec)%Init(iInit)%ParticleEmission &
!                = Species(iSpec)%Init(iInit)%PartDensity / Species(iSpec)%MacroParticleFactor * v_drift_line * A_ins
!            END IF
!            !---calculation of initial (macro)particle number
!            IF (Species(iSpec)%Init(iInit)%UseForInit) THEN
!              IF (Species(iSpec)%Init(iInit)%initialParticleNumber .GT. 0) THEN
!                CALL abort(&
!__STAMP__&
!                  ,'Either initialParticleNumber or PartDensity can be defined for selected parameters, not both!')
!              END IF
!              IF (TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'cuboid') THEN
!                Species(iSpec)%Init(iInit)%initialParticleNumber &
!                  = INT(Species(iSpec)%Init(iInit)%PartDensity / Species(iSpec)%MacroParticleFactor &
!                  * Species(iSpec)%Init(iInit)%CuboidHeightIC * A_ins)
!              ELSE !cylinder
!                Species(iSpec)%Init(iInit)%initialParticleNumber &
!                  = INT(Species(iSpec)%Init(iInit)%PartDensity / Species(iSpec)%MacroParticleFactor &
!                  * Species(iSpec)%Init(iInit)%CylinderHeightIC * A_ins)
!              END IF
!            END IF
!          ELSE
!            CALL abort(&
!__STAMP__&
!              ,'BaseVectors are parallel or zero!')
!          END IF
!        ELSE
!          CALL abort(&
!__STAMP__&
!          ,'Only const. or maxwell(_lpn) is supported as velocityDistr. for PartDensity without LD!')
!        END IF
!      ELSE IF (Species(iSpec)%Init(iInit)%VirtPreInsert) THEN
!        IF (Species(iSpec)%Init(iInit)%ParticleEmission .GT. 0.) THEN
!               CALL abort(&
!__STAMP__&
!          ,'Either ParticleEmission or PartDensity can be defined for selected emission parameters, not both!')
!        ELSE
!          SWRITE(*,*) "PartDensity is used for VPI of Species, Init ", iSpec, iInit !Value is calculated inside SetParticlePostion!
!        END IF
!      ELSE
!        CALL abort(&
!__STAMP__&
!        ,'PartDensity without LD is only supported for the SpaceIC cuboid(_vpi) or cylinder(_vpi)!')
!      END IF
!    END IF
!    !--- determine StartnumberOfInits (start loop index, e.g., for emission loops)
!    IF(iInit.EQ.0)THEN
!      !!!for new case: check if to be included!!!
!      IF((( (Species(iSpec)%Init(iInit)%initialParticleNumber.EQ.0)&
!        .AND.(Species(iSpec)%Init(iInit)%ParticleEmission.EQ.0.) )  &
!        .AND.(Species(iSpec)%Init(iInit)%PartDensity.EQ.0.) )       &
!        .AND.(Species(iSpec)%Init(iInit)%ConstantPressure.EQ.0.)    &
!        .AND.(Species(iSpec)%NumberOfInits.GT.0))       THEN 
!        Species(iSpec)%StartnumberOfInits = 1 ! only new style paramaters defined (Part-Species(i)-Init(iInit)-***)
!      ELSE
!        Species(iSpec)%StartnumberOfInits = 0 ! old style parameters has been defined for inits/emissions (Part-Species(i)-***)
!      END IF
!      SWRITE(*,*) "StartnumberOfInits of Species ", iSpec, " = ", Species(iSpec)%StartnumberOfInits
!    END IF ! iInit .EQ.0
!
  END DO ! iInit
END DO ! iSpec 
!
!! get information for IMD atom/ion charge determination and distribution
!IMDnSpecies         = GETINT('IMDnSpecies','1')
!IMDInputFile        = GETSTR('IMDInputFile','no file found')
!ALLOCATE(IMDSpeciesID(IMDnSpecies))
!ALLOCATE(IMDSpeciesCharge(IMDnSpecies))
!iIMDSpec=1
!DO iSpec = 1, nSpecies
!  WRITE(UNIT=hilf,FMT='(I2)') iSpec
!  IsIMDSpecies = GETLOGICAL('Part-Species'//TRIM(hilf)//'-IsIMDSpecies','.FALSE.')
!  IF(IsIMDSpecies)THEN
!    IMDSpeciesID(iIMDSpec)=iSpec
!    IMDSpeciesCharge(iIMDSpec)=NINT(Species(iSpec)%ChargeIC/1.60217653E-19)
!    iIMDSpec=iIMDSpec+1
!  END IF
!END DO
!
!
!! Which Lorentz boost method should be used?
!PartLorentzType = GETINT('Part-LorentzType','3')
!
! Read in boundary parameters
!dummy_int = CNTSTR('Part-nBounds')       ! check if Part-nBounds is present in .ini file
nPartBound = GETINT('Part-nBounds','1.')  ! get number of particle boundaries
IF ((nPartBound.LE.0).OR.(dummy_int.LT.0)) THEN
  CALL abort(&
__STAMP__&
  ,'ERROR: nPartBound .LE. 0:', nPartBound)
END IF
ALLOCATE(PartBound%SourceBoundName(1:nPartBound))
ALLOCATE(PartBound%TargetBoundCond(1:nPartBound))
ALLOCATE(PartBound%MomentumACC(1:nPartBound))
ALLOCATE(PartBound%WallTemp(1:nPartBound))
!ALLOCATE(PartBound%TransACC(1:nPartBound))
!ALLOCATE(PartBound%VibACC(1:nPartBound))
!ALLOCATE(PartBound%RotACC(1:nPartBound))
!ALLOCATE(PartBound%ElecACC(1:nPartBound))
!ALLOCATE(PartBound%Resample(1:nPartBound))
ALLOCATE(PartBound%WallVelo(1:3,1:nPartBound))
ALLOCATE(PartBound%AmbientCondition(1:nPartBound))
ALLOCATE(PartBound%AmbientConditionFix(1:nPartBound))
ALLOCATE(PartBound%AmbientTemp(1:nPartBound))
!ALLOCATE(PartBound%AmbientMeanPartMass(1:nPartBound))
!ALLOCATE(PartBound%AmbientBeta(1:nPartBound))
ALLOCATE(PartBound%AmbientVelo(1:3,1:nPartBound))
ALLOCATE(PartBound%AmbientDens(1:nPartBound))
!ALLOCATE(PartBound%AmbientDynamicVisc(1:nPartBound))
!ALLOCATE(PartBound%AmbientThermalCond(1:nPartBound))
ALLOCATE(PartBound%SolidState(1:nPartBound))
!ALLOCATE(PartBound%SolidCatalytic(1:nPartBound))
!ALLOCATE(PartBound%SolidSpec(1:nPartBound))
!ALLOCATE(PartBound%SolidPartDens(1:nPartBound))
!ALLOCATE(PartBound%SolidMassIC(1:nPartBound))
!ALLOCATE(PartBound%SolidAreaIncrease(1:nPartBound))
!ALLOCATE(PartBound%SolidCrystalIndx(1:nPartBound))
!ALLOCATE(PartBound%LiquidSpec(1:nPartBound))
!ALLOCATE(PartBound%ParamAntoine(1:3,1:nPartBound))
PartBound%SolidState(1:nPartBound)=.FALSE.
!PartBound%LiquidSpec(1:nPartBound)=0
!SolidSimFlag = .FALSE.
!LiquidSimFlag = .FALSE.
!
!ALLOCATE(PartBound%Adaptive(1:nPartBound))
!ALLOCATE(PartBound%AdaptiveType(1:nPartBound))
!ALLOCATE(PartBound%AdaptiveTemp(1:nPartBound))
!ALLOCATE(PartBound%AdaptivePressure(1:nPartBound))
!nAdaptiveBC = 0
!PartBound%Adaptive(:) = .FALSE.
!PartBound%AdaptiveType(:) = -1
!PartBound%AdaptiveTemp(:) = -1.
!PartBound%AdaptivePressure(:) = -1.
!
!ALLOCATE(PartBound%Voltage(1:nPartBound))
!ALLOCATE(PartBound%UseForQCrit(1:nPartBound))
!ALLOCATE(PartBound%Voltage_CollectCharges(1:nPartBound))
!PartBound%Voltage_CollectCharges(:)=0.
!ALLOCATE(PartBound%NbrOfSpeciesSwaps(1:nPartBound))
!!--determine MaxNbrOfSpeciesSwaps for correct allocation
!MaxNbrOfSpeciesSwaps=0
!DO iPartBound=1,nPartBound
!  WRITE(UNIT=hilf,FMT='(I2)') iPartBound
!  PartBound%NbrOfSpeciesSwaps(iPartBound)= GETINT('Part-Boundary'//TRIM(hilf)//'-NbrOfSpeciesSwaps','0')
!  MaxNbrOfSpeciesSwaps=max(PartBound%NbrOfSpeciesSwaps(iPartBound),MaxNbrOfSpeciesSwaps)
!END DO
!IF (MaxNbrOfSpeciesSwaps.gt.0) THEN
!  ALLOCATE(PartBound%ProbOfSpeciesSwaps(1:nPartBound))
!  ALLOCATE(PartBound%SpeciesSwaps(1:2,1:MaxNbrOfSpeciesSwaps,1:nPartBound))
!END IF
!--
!DO iPartBound=1,nPartBound
!  WRITE(UNIT=hilf,FMT='(I2)') iPartBound
!  tmpString = TRIM(GETSTR('Part-Boundary'//TRIM(hilf)//'-Condition','open'))
!  SELECT CASE (TRIM(tmpString))
!  CASE('open')
!     PartBound%TargetBoundCond(iPartBound) = PartBound%OpenBC          ! definitions see typesdef_pic
!     PartBound%AmbientCondition(iPartBound) = GETLOGICAL('Part-Boundary'//TRIM(hilf)//'-AmbientCondition','.FALSE.')
!     IF(PartBound%AmbientCondition(iPartBound)) THEN
!       PartBound%AmbientConditionFix(iPartBound) = GETLOGICAL('Part-Boundary'//TRIM(hilf)//'-AmbientConditionFix','.TRUE.')
!       PartBound%AmbientTemp(iPartBound) = GETREAL('Part-Boundary'//TRIM(hilf)//'-AmbientTemp','0')
!       PartBound%AmbientMeanPartMass(iPartBound) = GETREAL('Part-Boundary'//TRIM(hilf)//'-AmbientMeanPartMass','0')
!       PartBound%AmbientBeta(iPartBound) = &
!       SQRT(PartBound%AmbientMeanPartMass(iPartBound)/(2*BoltzmannConst*PartBound%AmbientTemp(iPartBound)))
!       PartBound%AmbientVelo(1:3,iPartBound) = GETREALARRAY('Part-Boundary'//TRIM(hilf)//'-AmbientVelo',3,'0. , 0. , 0.')
!       PartBound%AmbientDens(iPartBound) = GETREAL('Part-Boundary'//TRIM(hilf)//'-AmbientDens','0')
!       PartBound%AmbientDynamicVisc(iPartBound)=&
!           GETREAL('Part-Boundary'//TRIM(hilf)//'-AmbientDynamicVisc','1.72326582572253E-5') ! N2:T=288K
!       PartBound%AmbientThermalCond(iPartBound)=&
!           GETREAL('Part-Boundary'//TRIM(hilf)//'-AmbientThermalCond','2.42948500556027E-2') ! N2:T=288K
!     END IF
!     PartBound%Adaptive(iPartBound) = GETLOGICAL('Part-Boundary'//TRIM(hilf)//'-Adaptive','.FALSE.')
!     IF(PartBound%Adaptive(iPartBound)) THEN
!       nAdaptiveBC = nAdaptiveBC + 1
!       PartBound%AdaptiveType(iPartBound) = GETINT('Part-Boundary'//TRIM(hilf)//'-AdaptiveType','2')
!       PartBound%AdaptiveTemp(iPartBound) = GETREAL('Part-Boundary'//TRIM(hilf)//'-AdaptiveTemp','0.')
!       PartBound%AdaptivePressure(iPartBound) = GETREAL('Part-Boundary'//TRIM(hilf)//'-AdaptivePressure','0.')
!       IF (PartBound%AdaptiveTemp(iPartBound)*PartBound%AdaptivePressure(iPartBound).EQ.0.) THEN
!         CALL abort(&
!__STAMP__&
!,'Error during ParticleBoundary init: Part-Boundary'//TRIM(hilf)//'-AdaptiveTemp or -AdaptivePressure not defined')
!       END IF
!     END IF
!     PartBound%Voltage(iPartBound)         = GETREAL('Part-Boundary'//TRIM(hilf)//'-Voltage','0')
!  CASE('reflective')
!     PartBound%TargetBoundCond(iPartBound) = PartBound%ReflectiveBC
!     PartBound%MomentumACC(iPartBound)     = GETREAL('Part-Boundary'//TRIM(hilf)//'-MomentumACC','0')
!     PartBound%WallTemp(iPartBound)        = GETREAL('Part-Boundary'//TRIM(hilf)//'-WallTemp','0')
!     PartBound%TransACC(iPartBound)        = GETREAL('Part-Boundary'//TRIM(hilf)//'-TransACC','0')
!     PartBound%VibACC(iPartBound)          = GETREAL('Part-Boundary'//TRIM(hilf)//'-VibACC','0')
!     PartBound%RotACC(iPartBound)          = GETREAL('Part-Boundary'//TRIM(hilf)//'-RotACC','0')
!     PartBound%ElecACC(iPartBound)         = GETREAL('Part-Boundary'//TRIM(hilf)//'-ElecACC','0')
!     PartBound%Resample(iPartBound)        = GETLOGICAL('Part-Boundary'//TRIM(hilf)//'-Resample','.FALSE.')
!     PartBound%WallVelo(1:3,iPartBound)    = GETREALARRAY('Part-Boundary'//TRIM(hilf)//'-WallVelo',3,'0. , 0. , 0.')
!     PartBound%Voltage(iPartBound)         = GETREAL('Part-Boundary'//TRIM(hilf)//'-Voltage','0')
!     PartBound%SolidState(iPartBound)      = GETLOGICAL('Part-Boundary'//TRIM(hilf)//'-SolidState','.TRUE.')
!     PartBound%LiquidSpec(iPartBound)      = GETINT('Part-Boundary'//TRIM(hilf)//'-LiquidSpec','0')
!     IF(PartBound%SolidState(iPartBound))THEN
!       SolidSimFlag = .TRUE.
!       PartBound%SolidCatalytic(iPartBound)    = GETLOGICAL('Part-Boundary'//TRIM(hilf)//'-SolidCatalytic','.FALSE.')
!       PartBound%SolidSpec(iPartBound)         = GETINT('Part-Boundary'//TRIM(hilf)//'-SolidSpec','0')
!       PartBound%SolidPartDens(iPartBound)     = GETREAL('Part-Boundary'//TRIM(hilf)//'-SolidPartDens','1.0E+19')
!       PartBound%SolidMassIC(iPartBound)       = GETREAL('Part-Boundary'//TRIM(hilf)//'-SolidMassIC','3.2395E-25')
!       PartBound%SolidAreaIncrease(iPartBound) = GETREAL('Part-Boundary'//TRIM(hilf)//'-SolidAreaIncrease','1.')
!       PartBound%SolidCrystalIndx(iPartBound)  = GETINT('Part-Boundary'//TRIM(hilf)//'-SolidCrystalIndx','4')
!     END IF
!     IF (PartBound%LiquidSpec(iPartBound).GT.nSpecies) CALL abort(&
!__STAMP__&
!     ,'Particle Boundary Liquid Species not defined. Liquid Species: ',PartBound%LiquidSpec(iPartBound))
!     ! Parameters for evaporation pressure using Antoine Eq.
!     PartBound%ParamAntoine(1:3,iPartBound) = GETREALARRAY('Part-Boundary'//TRIM(hilf)//'-ParamAntoine',3,'0. , 0. , 0.')
!     IF ( (.NOT.PartBound%SolidState(iPartBound)) .AND. (ALMOSTZERO(PartBound%ParamAntoine(1,iPartBound))) &
!          .AND. (ALMOSTZERO(PartBound%ParamAntoine(2,iPartBound))) .AND. (ALMOSTZERO(PartBound%ParamAntoine(3,iPartBound))) ) THEN
!        CALL abort(&
!__STAMP__&
!       ,'Antoine Parameters not defined for Liquid Particle Boundary: ',iPartBound)
!     END IF
!     IF (.NOT.PartBound%SolidState(iPartBound)) LiquidSimFlag = .TRUE.
!     IF (PartBound%NbrOfSpeciesSwaps(iPartBound).gt.0) THEN  
!       !read Species to be changed at wall (in, out), out=0: delete
!       PartBound%ProbOfSpeciesSwaps(iPartBound)= GETREAL('Part-Boundary'//TRIM(hilf)//'-ProbOfSpeciesSwaps','1.')
!       DO iSwaps=1,PartBound%NbrOfSpeciesSwaps(iPartBound)
!         WRITE(UNIT=hilf2,FMT='(I2)') iSwaps
!         PartBound%SpeciesSwaps(1:2,iSwaps,iPartBound) = &
!             GETINTARRAY('Part-Boundary'//TRIM(hilf)//'-SpeciesSwaps'//TRIM(hilf2),2,'0. , 0.')
!       END DO
!     END IF
!  CASE('periodic')
!     PartBound%TargetBoundCond(iPartBound) = PartBound%PeriodicBC
!  CASE('simple_anode')
!     PartBound%TargetBoundCond(iPartBound) = PartBound%SimpleAnodeBC
!  CASE('simple_cathode')
!     PartBound%TargetBoundCond(iPartBound) = PartBound%SimpleCathodeBC
!  CASE('symmetric')
!     PartBound%TargetBoundCond(iPartBound) = PartBound%SymmetryBC
!     PartBound%WallVelo(1:3,iPartBound)    = (/0.,0.,0./)
!  CASE('analyze')
!     PartBound%TargetBoundCond(iPartBound) = PartBound%AnalyzeBC
!     IF (PartBound%NbrOfSpeciesSwaps(iPartBound).gt.0) THEN  
!       !read Species to be changed at wall (in, out), out=0: delete
!       PartBound%ProbOfSpeciesSwaps(iPartBound)= GETREAL('Part-Boundary'//TRIM(hilf)//'-ProbOfSpeciesSwaps','1.')
!       DO iSwaps=1,PartBound%NbrOfSpeciesSwaps(iPartBound)
!         WRITE(UNIT=hilf2,FMT='(I2)') iSwaps
!         PartBound%SpeciesSwaps(1:2,iSwaps,iPartBound) = &
!             GETINTARRAY('Part-Boundary'//TRIM(hilf)//'-SpeciesSwaps'//TRIM(hilf2),2,'0. , 0.')
!       END DO
!     END IF
!  CASE DEFAULT
!     SWRITE(*,*) ' Boundary does not exists: ', TRIM(tmpString)
!     CALL abort(&
!__STAMP__&
!         ,'Particle Boundary Condition does not exist')
!  END SELECT
!  PartBound%SourceBoundName(iPartBound) = TRIM(GETSTR('Part-Boundary'//TRIM(hilf)//'-SourceName'))
!  PartBound%UseForQCrit(iPartBound) = GETLOGICAL('Part-Boundary'//TRIM(hilf)//'-UseForQCrit','.TRUE.')
!  SWRITE(*,*)"PartBound",iPartBound,"is used for the Q-Criterion"
!END DO

!DEALLOCATE(PartBound%AmbientMeanPartMass)
!DEALLOCATE(PartBound%AmbientTemp)
!! Set mapping from field boundary to particle boundary index
!ALLOCATE(PartBound%MapToPartBC(1:nBCs))
!PartBound%MapToPartBC(:)=-10
!DO iPBC=1,nPartBound
!  DO iBC = 1, nBCs
!    IF (BoundaryType(iBC,1).EQ.0) THEN
!      PartBound%MapToPartBC(iBC) = -1 !there are no internal BCs in the mesh, they are just in the name list!
!      SWRITE(*,*)"... PartBound",iPBC,"is internal bound, no mapping needed"
!    END IF
!    IF (TRIM(BoundaryName(iBC)).EQ.TRIM(PartBound%SourceBoundName(iPBC))) THEN
!      PartBound%MapToPartBC(iBC) = iPBC !PartBound%TargetBoundCond(iPBC)
!      SWRITE(*,*)"... Mapped PartBound",iPBC,"on FieldBound",BoundaryType(iBC,1),",i.e.:",TRIM(BoundaryName(iBC))
!    END IF
!  END DO
!END DO
!! Errorhandler for PartBound-Types that could not be mapped to the 
!! FieldBound-Types.
!DO iBC = 1,nBCs
!  IF (PartBound%MapToPartBC(iBC).EQ.-10) THEN
!    CALL abort(&
!__STAMP__&
!    ,' PartBound%MapToPartBC for Boundary is not set. iBC: :',iBC)
!  END IF
!END DO

ALLOCATE(PEM%Element(1:PDM%maxParticleNumber), PEM%lastElement(1:PDM%maxParticleNumber), STAT=ALLOCSTAT) 
IF (ALLOCSTAT.NE.0) THEN
 CALL abort(&
__STAMP__&
  ,' Cannot allocate PEM arrays!')
END IF
!IF (useDSMC.OR.PartPressureCell) THEN
!  ALLOCATE(PEM%pStart(1:nElems)                         , &
!           PEM%pNumber(1:nElems)                        , &
!           PEM%pEnd(1:nElems)                           , &
!           PEM%pNext(1:PDM%maxParticleNumber)           , STAT=ALLOCSTAT) 
!           !PDM%nextUsedPosition(1:PDM%maxParticleNumber)
!  IF (ALLOCSTAT.NE.0) THEN
!    CALL abort(&
!__STAMP__&
!    , ' Cannot allocate DSMC PEM arrays!')
!  END IF
!END IF
!IF (useDSMC) THEN
!  ALLOCATE(PDM%PartInit(1:PDM%maxParticleNumber), STAT=ALLOCSTAT) 
!           !PDM%nextUsedPosition(1:PDM%maxParticleNumber)
!  IF (ALLOCSTAT.NE.0) THEN
!    CALL abort(&
!__STAMP__&
!    ,' Cannot allocate DSMC PEM arrays!')
!  END IF
!END IF

!--- initialize randomization (= random if one or more seeds are 0 or random is wanted)
!nrSeeds = GETINT('Part-NumberOfRandomSeeds','0')
!IF (nrSeeds.GT.0) THEN
!   ALLOCATE(seeds(1:nrSeeds))
!   DO iSeed = 1, nrSeeds
!      WRITE(UNIT=hilf,FMT='(I2)') iSeed
!      seeds(iSeed) = GETINT('Particles-RandomSeed'//TRIM(hilf),'0')
!   END DO
!END IF
!
!CALL RANDOM_SEED(Size = SeedSize)                       ! Check for number of needed Seeds
!TrueRandom = .FALSE.                             ! FALSE for defined random seed
!
!! to be stored in HDF5-state file!!!
!IF (nrSeeds.GT.0) THEN
!   IF (nrSeeds.NE.SeedSize) THEN
!      IF(printRandomSeeds)THEN
!        IPWRITE(UNIT_StdOut,*) 'Error: Number of seeds for RNG must be ',SeedSize
!        IPWRITE(UNIT_StdOut,*) 'Random RNG seeds are used'
!      END IF
!      TrueRandom = .TRUE.
!   END IF
!   DO iSeed = 1, nrSeeds
!      IF (Seeds(iSeed).EQ.0) THEN
!         IF(printRandomSeeds)THEN
!           IPWRITE(UNIT_StdOut,*) 'Error: ',SeedSize,' seeds for RNG must be defined'
!           IPWRITE(UNIT_StdOut,*) 'Random RNG seeds are used'
!         END IF
!         TrueRandom = .TRUE.
!      END IF
!   END DO
!ELSE
!   TrueRandom = .TRUE.
!END IF
!
!IF (TrueRandom) THEN
!   CALL RANDOM_SEED()
!ELSE
!#ifdef MPI
!   Seeds(1:SeedSize) = Seeds(1:SeedSize)+PartMPI%MyRank
!#endif
!   CALL RANDOM_SEED(PUT = Seeds(1:SeedSize))
!END IF
!
!ALLOCATE(iseeds(SeedSize))
!iseeds(:)=0
!CALL RANDOM_SEED(GET = iseeds(1:SeedSize))
!! to be stored in HDF5-state file!!!
!IF(printRandomSeeds)THEN
!  IPWRITE(UNIT_StdOut,*) 'Random seeds in PIC_init:'
!  DO iSeed = 1,SeedSize
!     IPWRITE(UNIT_StdOut,*) iseeds(iSeed)
!  END DO
!END IF
!DEALLOCATE(iseeds)
!DoZigguratSampling = GETLOGICAL('Particles-DoZigguratSampling','.FALSE.')
!DoPoissonRounding = GETLOGICAL('Particles-DoPoissonRounding','.FALSE.')
!DoTimeDepInflow   = GETLOGICAL('Particles-DoTimeDepInflow','.FALSE.')
!
!DO iSpec = 1, nSpecies
!  DO iInit = Species(iSpec)%StartnumberOfInits, Species(iSpec)%NumberOfInits
!    IF(Species(iSpec)%Init(iInit)%InflowRiseTime.GT.0.)THEN
!      IF(.NOT.DoPoissonRounding .AND. .NOT.DoTimeDepInflow)  CALL CollectiveStop(&
!__STAMP__, &
!' Linearly ramping of inflow-number-of-particles is only possible with PoissonRounding or DoTimeDepInflow!')
!    END IF      
!  END DO ! iInit = 0, Species(iSpec)%NumberOfInits
!END DO ! iSpec = 1, nSpecies


DelayTime = GETREAL('Part-DelayTime','0.')

!! init interpolation
!CALL InitializeInterpolation() ! not any more required ! has to be called earliear
!CALL InitPIC()
!! always, because you have to construct a halo_eps region around each bc element
!
#ifdef MPI
CALL MPI_BARRIER(PartMPI%COMM,IERROR)
#endif /*MPI*/
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)')' INIT FIBGM...' 
SafetyFactor  =GETREAL('Part-SafetyFactor','1.0')
halo_eps_velo =GETREAL('Particles-HaloEpsVelo','200')
!-- Finalizing InitializeVariables
CALL InitFIBGM()
!!CALL InitSFIBGM()
#ifdef MPI
CALL InitEmissionComm()
#endif /*MPI*/
#ifdef MPI
CALL MPI_BARRIER(PartMPI%COMM,IERROR)
#endif /*MPI*/

SWRITE(UNIT_StdOut,'(132("-"))')

!IF(enableParticleMerge) THEN
! CALL DefinePolyVec(vMPFMergePolyOrder) 
! CALL DefineSplitVec(vMPFMergeCellSplitOrder)
!END IF

END SUBROUTINE InitializeVariables

!----------------------------------------------------------------------------------------------------------------------------------!
! finalize particle variables
!----------------------------------------------------------------------------------------------------------------------------------!
SUBROUTINE FinalizeParticles() 
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Particle_Vars
!USE MOD_Particle_Mesh_Vars
USE MOD_Particle_Boundary_Vars
!!USE MOD_DSMC_Vars,                  ONLY: SampDSMC
!!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT VARIABLES 
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
!SDEALLOCATE(Pt_temp)
SDEALLOCATE(PartPosRef)
SDEALLOCATE(RandomVec)
SDEALLOCATE(PartState)
SDEALLOCATE(LastPartPos)
SDEALLOCATE(PartSpecies)
SDEALLOCATE(Pt)
SDEALLOCATE(PDM%ParticleInside)
SDEALLOCATE(PDM%nextFreePosition)
SDEALLOCATE(PDM%nextFreePosition)
SDEALLOCATE(PDM%dtFracPush)
SDEALLOCATE(PDM%IsNewPart)
!SDEALLOCATE(PartMPF)
!!SDEALLOCATE(Species%Init)
SDEALLOCATE(Species)
SDEALLOCATE(PartBound%SourceBoundName)
SDEALLOCATE(PartBound%TargetBoundCond)
SDEALLOCATE(PartBound%MomentumACC)
SDEALLOCATE(PartBound%WallTemp)
!SDEALLOCATE(PartBound%TransACC)
!SDEALLOCATE(PartBound%VibACC)
!SDEALLOCATE(PartBound%RotACC)
!SDEALLOCATE(PartBound%ElecACC)
!SDEALLOCATE(PartBound%Resample)
SDEALLOCATE(PartBound%WallVelo)
SDEALLOCATE(PartBound%AmbientCondition)
SDEALLOCATE(PartBound%AmbientConditionFix)
SDEALLOCATE(PartBound%AmbientTemp)
!SDEALLOCATE(PartBound%AmbientMeanPartMass)
!SDEALLOCATE(PartBound%AmbientBeta)
SDEALLOCATE(PartBound%AmbientVelo)
SDEALLOCATE(PartBound%AmbientDens)
!SDEALLOCATE(PartBound%AmbientDynamicVisc)
!SDEALLOCATE(PartBound%AmbientThermalCond)
!SDEALLOCATE(PartBound%Adaptive)
!SDEALLOCATE(PartBound%AdaptiveType)
!SDEALLOCATE(PartBound%AdaptiveTemp)
!SDEALLOCATE(PartBound%AdaptivePressure)
!SDEALLOCATE(PartBound%Voltage)
!SDEALLOCATE(PartBound%UseForQCrit)
!SDEALLOCATE(PartBound%Voltage_CollectCharges)
!SDEALLOCATE(PartBound%NbrOfSpeciesSwaps)
!SDEALLOCATE(PartBound%ProbOfSpeciesSwaps)
!SDEALLOCATE(PartBound%SpeciesSwaps)
!SDEALLOCATE(PartBound%MapToPartBC)
!SDEALLOCATE(PartBound%SolidState)
!SDEALLOCATE(PartBound%SolidCatalytic)
!SDEALLOCATE(PartBound%SolidSpec)
!SDEALLOCATE(PartBound%SolidPartDens)
!SDEALLOCATE(PartBound%SolidMassIC)
!SDEALLOCATE(PartBound%SolidAreaIncrease)
!SDEALLOCATE(PartBound%SolidCrystalIndx)
!SDEALLOCATE(PartBound%LiquidSpec)
!SDEALLOCATE(PartBound%ParamAntoine)
!SDEALLOCATE(PEM%Element)
!SDEALLOCATE(PEM%lastElement)
!SDEALLOCATE(PEM%pStart)
!SDEALLOCATE(PEM%pNumber)
!SDEALLOCATE(PEM%pEnd)
!SDEALLOCATE(PEM%pNext)
!SDEALLOCATE(seeds)
!SDEALLOCATE(RegionBounds)
!SDEALLOCATE(RegionElectronRef)
END SUBROUTINE FinalizeParticles

END MODULE MOD_ParticleInit
