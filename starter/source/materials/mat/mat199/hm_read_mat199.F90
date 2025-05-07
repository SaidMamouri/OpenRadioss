!Copyright>        OpenRadioss
!Copyright>        Copyright (C) 1986-2025 Altair Engineering Inc.
!Copyright>
!Copyright>        This program is free software: you can redistribute it and/or modify
!Copyright>        it under the terms of the GNU Affero General Public License as published by
!Copyright>        the Free Software Foundation, either version 3 of the License, or
!Copyright>        (at your option) any later version.
!Copyright>
!Copyright>        This program is distributed in the hope that it will be useful,
!Copyright>        but WITHOUT ANY WARRANTY; without even the implied warranty of
!Copyright>        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!Copyright>        GNU Affero General Public License for more details.
!Copyright>
!Copyright>        You should have received a copy of the GNU Affero General Public License
!Copyright>        along with this program.  If not, see <https://www.gnu.org/licenses/>.
!Copyright>
!Copyright>
!Copyright>        Commercial Alternative: Altair Radioss Software
!Copyright>
!Copyright>        As an alternative to this open-source version, Altair also offers Altair Radioss
!Copyright>        software under a commercial license.  Contact Altair to discuss further if the
!Copyright>        commercial version may interest you: https://www.altair.com/radioss/.
!copyright>        openradioss
! ======================================================================================================================

      !||====================================================================
      !||    hm_read_mat199_mod   ../starter/source/materials/mat/mat199/hm_read_mat199.F90
      !||--- called by ------------------------------------------------------
      !||    hm_read_mat          ../starter/source/materials/mat/hm_read_mat.F90
      !||====================================================================
      module hm_read_mat199_mod
      contains
  
      !||====================================================================
      !||    hm_read_mat199           ../starter/source/materials/mat/mat134/hm_read_mat199.F90
      !||--- called by ------------------------------------------------------
      !||    hm_read_mat              ../starter/source/materials/mat/hm_read_mat.F90
      !||--- calls      -----------------------------------------------------
      !||    hm_get_floatv            ../starter/source/devtools/hm_reader/hm_get_floatv.F
      !||    hm_option_is_encrypted   ../starter/source/devtools/hm_reader/hm_option_is_encrypted.F
      !||    init_mat_keyword         ../starter/source/materials/mat/init_mat_keyword.F
      !||--- uses       -----------------------------------------------------
      !||    elbuftag_mod             ../starter/share/modules1/elbuftag_mod.F
      !||    message_mod              ../starter/share/message_module/message_mod.F
      !||    submodel_mod             ../starter/share/modules1/submodel_mod.F
      !||====================================================================
      subroutine hm_read_mat199( mtag, matparam ,                              &
                   parmat   ,nuvar, unitab   ,lsubmodel,                       &
                   mat_id   ,titr     ,     iout       )
  ! ----------------------------------------------------------------------------------------------------------------------
  !                                                   modules
  ! ----------------------------------------------------------------------------------------------------------------------
      use elbuftag_mod
      use matparam_def_mod
      use unitab_mod
      use message_mod
      use submodel_mod
      use constant_mod , only : one ,two, zero,three,em20,six,pi,third
! --------------------------------------------------------------------------------------------------------------------
! --------------------------------------------------------------------------------------------------------------------
!                                                 implicit none
! --------------------------------------------------------------------------------------------------------------------
      implicit none
! ----------------------------------------------------------------------------------------------------------------------
!                                                   included files
! ----------------------------------------------------------------------------------------------------------------------
#include "my_real.inc"
      
!c-----------------------------------------------
!c   d u m m y   a r g u m e n t s
!c-----------------------------------------------
      integer, intent(in)                          :: mat_id
      integer, intent(in)                          :: iout
      integer, intent(out)                         :: nuvar 
      type (unit_type_),intent(in) ::unitab 
      type(submodel_data), dimension(nsubmod),intent(in) :: lsubmodel
      character(len=nchartitle) ,intent(in)             :: titr
      my_real, dimension(100)       ,intent(inout)   :: parmat  
      type(matparam_struct_) ,intent(inout) :: matparam
      type(mlaw_tag_), intent(inout)  :: mtag
!-----------------------------------------------
!   l o c a l   v a r i a b l e s
!-----------------------------------------------
      logical :: is_available,is_encrypted
      integer :: ilaw
      my_real :: rho0, young, bulk, shear 
      my_real :: e,nu,alpha,theta,eta3,beta3,eta6,beta6,sigy,a11,a12
      my_real :: phi_podogorski,phi_rosdehal,phi, inv_phir_p,inv_phir_r
!=======================================================================
      is_encrypted = .false.
      is_available = .false.
      ilaw  = 199
!--------------------------------------------------------
!
      call hm_option_is_encrypted(is_encrypted)
!
!--------------------------------------------------------
!     read input fields
!--------------------------------------------------------
      call hm_get_floatv('MAT_RHO'            ,rho0    ,is_available, lsubmodel, unitab)
      !line2
      call hm_get_floatv('MAT_E'        ,e   ,is_available, lsubmodel, unitab)
      call hm_get_floatv('MAT_NU'       ,nu  ,is_available, lsubmodel, unitab)
      call hm_get_floatv('MAT_THETA'    ,theta  ,is_available, lsubmodel, unitab)
      call hm_get_floatv('MAT_ALPHA'    ,alpha  ,is_available, lsubmodel, unitab)
      call hm_get_floatv('MAT_ETA3'     ,eta3  ,is_available, lsubmodel, unitab)
      call hm_get_floatv('MAT_BETA3'    ,beta3 ,is_available, lsubmodel, unitab)
      call hm_get_floatv('MAT_ETA6'     ,eta6  ,is_available, lsubmodel, unitab)
      call hm_get_floatv('MAT_BETA6'    ,beta6 ,is_available, lsubmodel, unitab)
      call hm_get_floatv('MAT_SIGY'     ,sigy ,is_available, lsubmodel, unitab)
!-------------------------------------
      young  = e
      shear  = young/(two * (one + nu))
      bulk   = young/(three * (one - two * nu))
      a11    =  e*(one-nu)/(one-two*nu)/(one + nu)
      a12    = e*nu/(one-two*nu)/(one + nu)
      !
      inv_phir_p = cos(third* (pi * beta3 - acos(eta3)))
      inv_phir_r = cos((one / six) * (pi * beta6 - acos(eta6))) 
      if(inv_phir_p /=0) then
        inv_phir_p = one / inv_phir_p
      else
        inv_phir_p = em20
      endif
      if(inv_phir_r /=0) then
        inv_phir_r = one / inv_phir_r
      else
        inv_phir_r = em20
      endif
      phi_podogorski = inv_phir_p*cos(third* (pi * beta3 - acos(eta3 * cos(three* theta))))
      phi_rosdehal = inv_phir_r*cos((one/ six) * (pi * beta6 - acos(eta6 * cos(six* theta)))) 
      phi = alpha*phi_podogorski + (one - alpha)*phi_rosdehal
!-------------------------------------
      nuvar = 6
!-------------------------------------
      
      matparam%niparam = 0
      matparam%nuparam = 12
      matparam%nfunc   = 0
      matparam%ntable  = 0
!          
      allocate (matparam%uparam(matparam%nuparam))
      allocate (matparam%iparam(matparam%niparam))
      allocate (matparam%table(matparam%ntable))
!     
      matparam%uparam(1)  = e
      matparam%uparam(2)  = nu    
      matparam%uparam(3)  = theta
      matparam%uparam(4)  = alpha
      matparam%uparam(5)  = eta3
      matparam%uparam(6)  = beta3
      matparam%uparam(7)  = eta6
      matparam%uparam(8)  = beta6
      matparam%uparam(9)  = sigy
      matparam%uparam(10) = a11
      matparam%uparam(11) = a12
      matparam%uparam(12) = phi
      !
       !< Real material parameters
      matparam%rho       = rho0
      matparam%rho0      = rho0
      matparam%young      = young
      matparam%nu         = nu
      matparam%shear      = shear
      matparam%bulk       = bulk
      !
      parmat(1) = bulk
      parmat(2) = young
      parmat(3) = nu 
      parmat(16) = 1
      parmat(17) = (one-two*nu)/(one-nu)  !   2G / (bulk + G*4/3)
      !
      !mtag%g_epsd = 1 ! not used 
      mtag%g_pla   = 1   
      !MTAG%G_DMG   = 1
       !
      !MTAG%L_EPSD  = 1  
      mtag%l_pla   = 1  
      !MTAG%L_DMG   = 1
!-------------------------------------------------
              !< Properties compatibility  
          call init_mat_keyword(matparam,"SOLID_ISOTROPIC") 
          !< Properties compatibility  
          call init_mat_keyword(matparam ,"INCREMENTAL" )
          call init_mat_keyword(matparam ,"HOOK")
          call init_mat_keyword(matparam ,"ISOTROPIC")  

          call init_mat_keyword(matparam ,"SPH")  
!-------------------------------------------------
      write(iout,1050) trim(titr),mat_id,199
      write(iout,1000)
      if (is_encrypted) then
        write(iout,'(5x,a,//)')'CONFIDENTIAL DATA'
      else
        write(iout,1060) rho0
        write(iout,1100) e,nu,alpha,theta,eta3,beta3,eta6,beta6,sigy
      endif       
!     
!-----------
      return
!-----------
 1000 format(                                                                &
     5x,a,/,                                                                 & 
     5x,40h  GINIRALIZED MATERIAL                  ,/,                       & 
     5x,40h  -----------------------------------   ,//)           
 1050 format(/                                                               &
      5x,a,/,                                                                &
      5x,'MATERIAL NUMBER . . . . . . . . . . . . .=',i10/,                  &
      5x,'MATERIAL LAW. . . . . . . . . . . . . . .=',i10/)       
 1060 format(                                                                &
      5x,'INITIAL DENSITY . . . . . . . . . . . . .=',1pg20.13/)  
 1100 format(                                                                & 
      5x,'YOUNGS MODULUS . . . . . . . . . . . . . .  . . . ..=',1PG20.13/,  &
      5x,'POISSONS RATIO . . . . . . . . . . . . . . . . . . .=',1PG20.13,   &  
      5x,'THETA. . . . . . . . . . . . . . . . . . . . . . . .=',1PG20.13/,  &
      5x,'ALPHA . . .. . . . . . . . . . . . . . . . . . . . .=',1PG20.13/,  &
      5x,'ETA3  .  . . . . . . . . . . . . . . . . . . . . .. =',1PG20.13/,  &
      5x,'BETA3 . . . . . . . . .. . . . . . . . . . . . . . .=',1PG20.13/,  &
      5x,'ETA6  .  . . . . . . . . . . . . . . . . . . . . .. =',1PG20.13/,  &
      5x,'BETA6 . . . . . . . . .. . . . . . . . . . . . . . .=',1PG20.13/,  &     
      5x,'SIGY . . . . . . . . .. . . . . . . . . . . . . . ..=',1PG20.13/)
!-----------------
    end subroutine hm_read_mat199
end module hm_read_mat199_mod             
