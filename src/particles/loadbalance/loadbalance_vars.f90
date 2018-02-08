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
!===================================================================================================================================
! Variables needed for the evaluation of the record points
!===================================================================================================================================
MODULE MOD_LoadBalance_Vars
! MODULES
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
!LOGICAL                             :: DoLoadBalance                              ! DoLoadBalance
!LOGICAL                             :: InitLoadBalanceIsDone                      ! switch for checking
! time measurement
!REAL,ALLOCATABLE                    :: tTotal(:)                                  ! time measurement over whole dt_analyze 
!REAL,ALLOCATABLE                    :: tCurrent(:)                                ! time measurement over one step
!REAL,ALLOCATABLE                    :: LoadSum(:)                                 ! sum of load per step over whole dt_analyze 
!REAL(KIND=8)                        :: nTotalParts                                ! number of particles in time of tTotal
!INTEGER                             :: nLoadIter                                  ! number of load iter 
!!INTEGER                             :: nCurrentParts                              ! number of current particles
!INTEGER                             :: nLoadBalance                               ! number of load balances
!INTEGER                             :: nLoadBalanceSteps                          ! number of performed  load balances steps
!REAL,ALLOCATABLE                    :: LoadDistri(:)                              ! Weighted load distribution of all procs
!INTEGER,ALLOCATABLE                 :: PartDistri(:)                              ! Part distribution of all procs
!INTEGER                             :: PartWeightMethod                           ! method to compute the particle weight
!INTEGER                             :: WeightAverageMethod                        ! method to average the particle weight
!                                                                                  ! (1: iter, 2: dt_Analyze)
!                                                                                  ! nSkipAnalyze is greater than 1
!-----------------------------------------------------------------------------------------------------------------------------------
! particle load balancing
!-----------------------------------------------------------------------------------------------------------------------------------
!INTEGER                             :: nSkipAnalyze                               ! Skip Analyze-Dt
!REAL                                :: ParticleMPIWeight
!REAL                                :: DeviationThreshold                         ! threshold for load-balancing
!LOGICAL                             :: writePartitionInfo                         ! write partitioninfo file
!REAL                                :: WeightSum                                  ! global sum of all weights
!REAL                                :: targetWeight                               ! optimal weight for each proc
!-----------------------------------------------------------------------------------------------------------------------------------
! Element Local measurement
!-----------------------------------------------------------------------------------------------------------------------------------
!REAL                                :: tCartMesh                                  ! time for CartMesh deposition
!REAL                                :: tTracking                                  ! time for relocation of particles
!REAL,ALLOCATABLE                    :: ElemTime(:)
!REAL,ALLOCATABLE                    :: ElemGlobalTime(:)
INTEGER(KIND=8),ALLOCATABLE         :: nPartsPerElem(:)
!INTEGER(KIND=8),ALLOCATABLE         :: nDeposPerElem(:)
!INTEGER(KIND=8),ALLOCATABLE         :: nTracksPerElem(:)


END MODULE MOD_LoadBalance_Vars