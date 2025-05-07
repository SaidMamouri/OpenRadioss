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
!===============================================================================

      !||====================================================================
      !||    sigeps199s_mod   ../engine/source/materials/mat/mat134/sigeps199s.F90
      !||--- called by ------------------------------------------------------
      !||    mulaw            ../engine/source/materials/mat_share/mulaw.F90
      !||====================================================================
      module sigeps199s_mod
      contains
      !||====================================================================
      !||    sigeps199s         ../engine/source/materials/mat/mat199/sigeps199s.F90
      !||--- called by ------------------------------------------------------
      !||    mulaw              ../engine/source/materials/mat_share/mulaw.F90
      !||--- uses       -----------------------------------------------------
      !||    constant_mod       ../common_source/modules/constant_mod.F
      !||    matparam_def_mod   ../common_source/modules/mat_elem/matparam_def_mod.F90
      !||====================================================================
       subroutine sigeps199s(mat_param  ,                                       &
           nel      ,timestep ,                                                 &
           depsxx   ,depsyy   ,depszz   ,depsxy   ,depsyz   ,depszx   ,         &
           sigoxx   ,sigoyy   ,sigozz   ,sigoxy   ,sigoyz   ,sigozx   ,         &
           signxx   ,signyy   ,signzz   ,signxy   ,signyz   ,signzx   ,         &
           yld      ,pla      ,dpla     , etse     ,soundsp  ,off      )
!
! =================================================================================
! \brief orthotropic hill material with plastic strain rate dependancy for solids

! =================================================================================
!   m o d u l e s
!-----------------------------------------------
      use matparam_def_mod
      use constant_mod ,only : pi,zero,one,half,third,two,three,four,em20,six
      use constant_mod ,only : four_over_3,four_over_5,em10
! ---------------------------------------------------------------------------------
          implicit none
! ---------------------------------------------------------------------------------
!     included files
! ---------------------------------------------------------------------------------

#include "my_real.inc"

!-----------------------------------------------
!   d u m m y   a r g u m e n t s
!-----------------------------------------------
      integer ,intent(in) :: nel                           !< element group size
     !! integer ,intent(in) :: nuvar                         !< number of state variables
      my_real ,intent(in) :: timestep                      !< time step
      ! my_real ,dimension(nel)     ,intent(in)    :: rho0  !< reference density
      !!my_real ,dimension(nel)     ,intent(in)    :: rho   !< density 
      my_real ,dimension(nel)     ,intent(in)    :: depsxx !<  strain increment component in direction xx 
      my_real ,dimension(nel)     ,intent(in)    :: depsyy !< strain increment component in direction yy
      my_real ,dimension(nel)     ,intent(in)    :: depszz !< strain increment component in direction zz
      my_real ,dimension(nel)     ,intent(in)    :: depsxy !< strain increment component  in xy direction
      my_real ,dimension(nel)     ,intent(in)    :: depsyz !<   
      my_real ,dimension(nel)     ,intent(in)    :: depszx !< strain rate component  in zx direction
      my_real ,dimension(nel)     ,intent(in)    :: sigoxx !< output stress component
      my_real ,dimension(nel)     ,intent(in)    :: sigoyy !< output stress component
      my_real ,dimension(nel)     ,intent(in)    :: sigozz !< output stress component
      my_real ,dimension(nel)     ,intent(in)    :: sigoxy !< output stress component
      my_real ,dimension(nel)     ,intent(in)    :: sigoyz !< output stress component
      my_real ,dimension(nel)     ,intent(in)    :: sigozx !< output stress component
      my_real ,dimension(nel)     ,intent(out)   :: signxx !< output stress component
      my_real ,dimension(nel)     ,intent(out)   :: signyy !< output stress component
      my_real ,dimension(nel)     ,intent(out)   :: signzz !< output stress component
      my_real ,dimension(nel)     ,intent(out)   :: signxy !< output stress component
      my_real ,dimension(nel)     ,intent(out)   :: signyz !< output stress component
      my_real ,dimension(nel)     ,intent(out)   :: signzx !< output stress component
      my_real ,dimension(nel)     ,intent(inout) :: yld     !< yield function
      my_real ,dimension(nel)     ,intent(inout) :: pla     !< plastic strain
      my_real ,dimension(nel)     ,intent(out)   :: dpla    !< incremental plastic strain 
      my_real ,dimension(nel)     ,intent(out)   :: etse    !< 
      !
      !
      my_real ,dimension(nel)     ,intent(inout) :: off    !< element activation coefficient
      my_real ,dimension(nel)     ,intent(out)   :: soundsp!< sound speed
     !! my_real ,dimension(nel,nuvar)   ,intent(inout) :: uvar      !< state variables
      type (matparam_struct_)         ,intent(in)    :: mat_param !< material parameter structure
      target :: mat_param
!-----------------------------------------------
!   l o c a l   v a r i a b l e s
!-----------------------------------------------
      integer :: i,ndex,indx(nel),j
      my_real ,dimension(nel) :: p,sigeq,sxx,syy,szz
      my_real ::   a11,a12,alpha,eta3,beta3,eta6,beta6,sigy
      my_real ::   dtime,nu,young,shear,bulk,rho0,g3,theta,i2
      my_real ::  phi_podogorski,phi_rosdehal,scale,phi
!===============================================================================    
      dtime  = max(timestep, em20)
      young  = mat_param%young
      shear  = mat_param%shear
      bulk   = mat_param%bulk
      nu     = mat_param%nu
      rho0   = mat_param%rho0
      g3     = three*shear
!
      theta = mat_param%uparam(3)  
      alpha = mat_param%uparam(4)  
      eta3  = mat_param%uparam(5)  
      beta3 = mat_param%uparam(6)  
      eta6  = mat_param%uparam(7)  
      beta6 = mat_param%uparam(8)  
      sigy  = mat_param%uparam(9)  
      a11   = mat_param%uparam(10) 
      a12   = mat_param%uparam(11) 
      phi   = mat_param%uparam(12)
     !
      soundsp(1:nel) = sqrt((bulk + four_over_3*shear) / rho0)     ! sound-speed
      etse(1:nel)  = one
      yld(1:nel)  = sigy
!---------------------------------------------------------------------
      !  estimated elastic stress
       ndex = 0
       do i=1,nel
        ! elastic stress 
        signxx(i)  = sigoxx(i) + a11*depsxx(i) + a12*(depsyy(i) + depszz(i))
        signyy(i)  = sigoyy(i) + a11*depsyy(i) + a12*(depsxx(i) + depszz(i)) 
        signzz(i)  = sigozz(i) + a11*depszz(i) + a12*(depsxx(i) + depsyy(i)) 
        signxy(i)  = sigoxy(i) + shear*depsxy(i)
        signyz(i)  = sigoyz(i) + shear*depsyz(i)
        signzx(i)  = sigozx(i) + shear*depszx(i)
        ! deviatoric stress
         p(i) = third*(signxx(i) + signyy(i) + signzz(i))
         sxx(i) = signxx(i) - p(i) 
         syy(i) = signyy(i) - p(i)
         szz(i) = signzz(i) - p(i)
        ! computing equivalent stress 
         i2 = half*(sxx(i)**2 + syy(i)**2 + szz(i)**2) +       &
                              (signxy(i)**2 + signyz(i)**2 + signzx(i)**2)
         !! phi_podogorski = inv_phir_p*cos(third* (pi * beta3 - acos(eta3 * cos(three* theta))))
         !! phi_rosdehal = inv_phir_r*cos((one/ six) * (pi * beta6 - acos(eta6 * cos(six* theta)))) 
         sigeq(i) = sqrt(three * i2) * phi
         if(sigeq(i) >= sigy) then
           ndex = ndex + 1
           indx(ndex) = i
         endif
      enddo 
      !
      !  plasticity
      if(ndex == 0) return
      !  plasticity
      do j=1,ndex
         i = indx(j)
         scale  = sigy/sigeq(i) 
         signxx(i) = sxx(i)*scale + p(i)
         signyy(i) = syy(i)*scale + p(i)
         signzz(i) = szz(i)*scale + p(i)
         signxy(i) = signxy(i)*scale
         signyz(i) = signyz(i)*scale
         signzx(i) = signzx(i)*scale
         ! equivalent plastic strain
         dpla(i) = (one - scale)*sigeq(i)/max(g3,em20)
         pla(i) = pla(i) + dpla(i)
      enddo  
      do i=1,nel
        signxx(i) = signxx(i)*off(i)
        signyy(i) = signyy(i)*off(i)
        signzz(i) = signzz(i)*off(i)
        signxy(i) = signxy(i)*off(i)
        signyz(i) = signyz(i)*off(i)
        signzx(i) = signzx(i)*off(i)
      enddo  
!-----------
      return
      end
!-----------
      end module sigeps199s_mod
