Copyright>        OpenRadioss
Copyright>        Copyright (C) 1986-2025 Altair Engineering Inc.
Copyright>
Copyright>        This program is free software: you can redistribute it and/or modify
Copyright>        it under the terms of the GNU Affero General Public License as published by
Copyright>        the Free Software Foundation, either version 3 of the License, or
Copyright>        (at your option) any later version.
Copyright>
Copyright>        This program is distributed in the hope that it will be useful,
Copyright>        but WITHOUT ANY WARRANTY; without even the implied warranty of
Copyright>        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
Copyright>        GNU Affero General Public License for more details.
Copyright>
Copyright>        You should have received a copy of the GNU Affero General Public License
Copyright>        along with this program.  If not, see <https://www.gnu.org/licenses/>.
Copyright>
Copyright>
Copyright>        Commercial Alternative: Altair Radioss Software
Copyright>
Copyright>        As an alternative to this open-source version, Altair also offers Altair Radioss
Copyright>        software under a commercial license.  Contact Altair to discuss further if the
Copyright>        commercial version may interest you: https://www.altair.com/radioss/.
      !||====================================================================
      !||    arezon                 ../engine/source/ale/arezon.F
      !||--- called by ------------------------------------------------------
      !||    alethe                 ../engine/source/ale/alethe.F
      !||--- calls      -----------------------------------------------------
      !||    arezo2                 ../engine/source/ale/ale2d/arezo2.F
      !||    arezo2_sig             ../engine/source/ale/ale2d/arezo2.F
      !||    arezo3                 ../engine/source/ale/ale3d/arezo3.F
      !||    arezo3_sig             ../engine/source/ale/ale3d/arezo3.F
      !||    brezo2                 ../engine/source/ale/ale2d/brezo2.F
      !||    brezo2_sig             ../engine/source/ale/ale2d/brezo2.F
      !||    initbuf                ../engine/share/resol/initbuf.F
      !||    my_barrier             ../engine/source/system/machine.F
      !||    nvar                   ../engine/source/input/nvar.F
      !||    spmd_e1vois            ../engine/source/mpi/fluid/spmd_cfd.F
      !||--- uses       -----------------------------------------------------
      !||    ale_connectivity_mod   ../common_source/modules/ale/ale_connectivity_mod.F
      !||    elbufdef_mod           ../common_source/modules/mat_elem/elbufdef_mod.F90
      !||    initbuf_mod            ../engine/share/resol/initbuf.F
      !||    multimat_param_mod     ../common_source/modules/multimat_param_mod.F90
      !||====================================================================
      SUBROUTINE AREZON(
     1  IPARG   ,ELBUF_STR,FLUX    ,PHI   ,ALE_CONNECT   ,
     2  NDIM    ,NVAR     ,NV      ,ITASK ,NERCVOIS,
     3  NESDVOIS,LERCVOIS ,LESDVOIS,LENCOM,BHOLE   ,
     4  ITRIMAT ,IXS)
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      USE INITBUF_MOD
      USE ELBUFDEF_MOD            
      USE ALE_CONNECTIVITY_MOD
      USE MULTIMAT_PARAM_MOD , ONLY : M51_N0PHAS, M51_NVPHAS
C-----------------------------------------------
C   I m p l i c i t   T y p e s
C-----------------------------------------------
#include      "implicit_f.inc"
C-----------------------------------------------
C   C o m m o n   B l o c k s
C-----------------------------------------------
#include      "com01_c.inc"
#include      "com04_c.inc"
#include      "vect01_c.inc"
#include      "param_c.inc"
#include      "task_c.inc"
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      my_real FLUX(*), PHI(*)
      INTEGER NDIM, NVAR, ITASK, LENCOM,ITRIMAT,NV,
     .        IPARG(NPARG,NGROUP), 
     .        NERCVOIS(*),NESDVOIS(*),LERCVOIS(*),LESDVOIS(*),
     .        BHOLE(*),IXS(NIXS,NUMELS)
      TYPE (ELBUF_STRUCT_), DIMENSION(NGROUP), TARGET :: ELBUF_STR
      TYPE(t_ale_connectivity), INTENT(IN) :: ALE_CONNECT
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: NG, IRS, IRE, I, J, K, NM, NMN, NFX, ADD, ADD0,IDX,INDX,NEL
      my_real, DIMENSION(:), POINTER :: VAR,SIG,VOL, TAG22,TEMP
C-----------------------------------------------
C     S o u r c e   L i n e s
C-----------------------------------------------

C=======================================================================
C PHI DEFINED DEPENDING ON NVAR VALUE
C=======================================================================
      IDX=NV
      NMN=MAX(1,NMULT)
      DO NM=1,NMN
C---------------------
       CALL MY_BARRIER
C---------------------
       DO NG=ITASK+1,NGROUP,NTHREAD
C     ALE ON / OFF
        IF (IPARG(76, NG) == 1) CYCLE ! --> OFF
         CALL INITBUF(IPARG    ,NG      ,                  
     2     MTN     ,LLT     ,NFT     ,IAD     ,ITY     ,   
     3     NPT     ,JALE    ,ISMSTR  ,JEUL    ,JTUR    ,   
     4     JTHE    ,JLAG    ,JMULT   ,JHBE    ,JIVF    ,   
     5     NVAUX   ,JPOR    ,JCVT    ,JCLOSE  ,JPLASOL ,   
     6     IREP    ,IINT    ,IGTYP   ,ISRAT   ,ISROT   ,   
     7     ICSEN   ,ISORTH  ,ISORTHG ,IFAILURE,JSMS    )
        !-----------------------------!
        ! UNPLUG CONDITIONS           !
        !-----------------------------!
        IF (NVAR == 11 .AND. MTN /= 81) CYCLE
        IF (JALE+JEUL == 0)             CYCLE
        IF (IPARG(8,NG) == 1)           CYCLE
        IF (MAX(1,JMULT) < NM)         CYCLE
        ! pressurer rezoning for outlets (continuity)
        IF( JMULT /= 0) MTN =IPARG(24+NM,NG)
        IF (NVAR == 10 .AND.(MTN == 37)) CYCLE
        IF (NVAR == 12 .AND. ELBUF_STR(NG)%GBUF%G_TEMP == 0) CYCLE
        !-----------------------------!
        IRS=IPARG(15,NG)
        IRE=IPARG(16,NG)
        LFT=1
        NEL=LLT

        !----------------------------!
        !        N V A R = 2         !
        !            (SIGMA)         !
        !----------------------------!
        IF (NVAR == 2 .AND. IRS == 1) THEN
          IF (ITRIMAT>0.AND.MTN == 51) THEN
            ADD    = M51_N0PHAS + (ITRIMAT-1)*M51_NVPHAS + IDX          ! ADD+1 => SIG(1)
            ADD    = ADD *LLT
            DO I=LFT,LLT
              J      = I+NFT
              PHI(J) = ELBUF_STR(NG)%BUFLY(1)%MAT(1,1,1)%VAR(ADD+I)
            END DO
          ELSE
            DO I=LFT,LLT !other material laws or mtn=51
              J = I+NFT
              K = (IDX-1)*NEL + I
              PHI(J) = ELBUF_STR(NG)%BUFLY(NM)%LBUF(1,1,1)%SIG(K)
            ENDDO
          ENDIF
        !----------------------------!
        !        N V A R = 10        !
        !              (PLAS)        !
        !                (TB) explo. !
        !                (RK) turb.  !
        !----------------------------!
        ELSEIF (NVAR == 10 .AND. IRE == 1) THEN
           IF (MTN == 41) CYCLE
          IF (MTN==51 .AND. ITRIMAT /=4) THEN  !submat 4 is JWL
            IF(ITRIMAT==0)CYCLE
            ADD0= M51_N0PHAS + (ITRIMAT-1)*M51_NVPHAS
            ADD = ADD0 + 15  ! PLAS
            K   = LLT*(ADD-1)
            VAR  => ELBUF_STR(NG)%BUFLY(1)%MAT(1,1,1)%VAR(K+1:K+LLT*NDIM)
          ELSEIF (MTN == 5 .OR. MTN ==97 .OR. MTN==105 .OR. ITRIMAT == 4) THEN  ! DETONATION TIME
            VAR => ELBUF_STR(NG)%BUFLY(NM)%LBUF(1,1,1)%TB(1:LLT)
          ELSEIF (MTN == 6) THEN
            VAR => ELBUF_STR(NG)%BUFLY(NM)%LBUF(1,1,1)%RK(1:LLT)
          ELSEIF (MTN >= 28 .AND. MTN /= 67 .AND. MTN /= 49) THEN
            VAR => ELBUF_STR(NG)%BUFLY(NM)%LBUF(1,1,1)%PLA(1:LLT)
          ELSE
            VAR => ELBUF_STR(NG)%BUFLY(NM)%LBUF(1,1,1)%PLA(1:LLT)
          ENDIF
          DO I=LFT,LLT
            J = I+NFT
            PHI(J) = VAR(I)
          ENDDO
        !----------------------------!
        !        N V A R = 11        !
        !             (EPS_P_VOL)    !
        !----------------------------!
        ELSEIF (NVAR == 11 .AND. IRE == 1 .AND. MTN == 81) THEN
          VAR => ELBUF_STR(NG)%BUFLY(NM)%MAT(1,1,1)%VAR(LLT+1:LLT*2)
          DO I=LFT,LLT
            J = I+NFT
            PHI(J) = VAR(I)
          ENDDO
        !----------------------------!
        !        N V A R = 12        !
        !             (TEMP)         !
        !----------------------------!
        ELSEIF (NVAR == 12) THEN
          VAR => ELBUF_STR(NG)%BUFLY(NM)%LBUF(1,1,1)%TEMP(1:LLT)
          DO I=LFT,LLT
            J = I+NFT
            PHI(J) = VAR(I)
          ENDDO
        !----------------------------!
        !        DEFAULT             !
        !----------------------------!
        ELSE
          DO I=LFT,LLT
            J=I+NFT
            PHI(J)=ZERO
          ENDDO
        ENDIF !(NVAR == 2 .AND. IRS == 1)
      ENDDO  ! NG
C---------------------
       CALL MY_BARRIER
C---------------------
       !------------------------------------!
       ! SPMD EXCHANGES                     !
       !------------------------------------!
       IF (NSPMD > 1)THEN
!$OMP SINGLE
         CALL SPMD_E1VOIS(PHI,NERCVOIS,NESDVOIS,LERCVOIS,LESDVOIS,LENCOM  )
!$OMP  END SINGLE
       END IF


C=======================================================================
C REZONING DEPENDING ON NVAR VALUE
C=======================================================================

       DO NG=ITASK+1,NGROUP,NTHREAD
C     ALE ON / OFF
        IF (IPARG(76, NG) == 1) CYCLE ! --> OFF
        CALL INITBUF(IPARG    ,NG      ,
     2     MTN     ,LLT     ,NFT     ,IAD     ,ITY     ,
     3     NPT     ,JALE    ,ISMSTR  ,JEUL    ,JTUR    ,
     4     JTHE    ,JLAG    ,JMULT   ,JHBE    ,JIVF    ,
     5     NVAUX   ,JPOR    ,JCVT    ,JCLOSE  ,JPLASOL ,
     6     IREP    ,IINT    ,IGTYP   ,ISRAT   ,ISROT   ,
     7     ICSEN   ,ISORTH  ,ISORTHG ,IFAILURE,JSMS    )
        !-----------------------------!
        ! UNPLUG CONDITIONS           !
        !-----------------------------!
        IF (NVAR == 11 .AND. MTN /= 81) CYCLE
        IF (MAX(1,JMULT) < NM)         CYCLE
        IF (JALE+JEUL == 0)             CYCLE
        IF (IPARG(8,NG) == 1)           CYCLE
        IF (ITRIMAT /= 0.AND.MTN /= 51) CYCLE
        IF (JMULT /= 0) MTN = IPARG(24+NM,NG)
        IF (NVAR == 10.AND.(MTN == 37)) CYCLE
        IF (NVAR == 10.AND.MTN==51.AND.ITRIMAT==0) CYCLE
        IF (NVAR == 12 .AND. ELBUF_STR(NG)%GBUF%G_TEMP == 0) CYCLE
        !-----------------------------!
        IRS=IPARG(15,NG)
        IRE=IPARG(16,NG)
        NEL=LLT
        !----------------------------!
        !        N V A R = 2         !
        !            (SIGMA)         !
        !----------------------------!
        IF (NVAR == 2 .AND. IRS == 1) THEN
          INDX = IDX
          IF (ITRIMAT>0) THEN
            ADD    = M51_N0PHAS + (ITRIMAT-1)*M51_NVPHAS + IDX          ! ADD+1 => SIG[IDX](i=1)
            ADD    = ADD *LLT
            SIG    => ELBUF_STR(NG)%BUFLY(1)%MAT(1,1,1)%VAR(ADD+1:ADD+LLT*NDIM)
            ADD    = M51_N0PHAS + (ITRIMAT-1)*M51_NVPHAS + 10          ! ADD+1 => VOL(1)
            ADD    = ADD *LLT
            VOL => ELBUF_STR(NG)%BUFLY(1)%MAT(1,1,1)%VAR(ADD+1:ADD+LLT)
            INDX   = 1
          ELSE
            SIG => ELBUF_STR(NG)%BUFLY(NM)%LBUF(1,1,1)%SIG(1:LLT*NDIM)
            VOL => ELBUF_STR(NG)%BUFLY(NM)%LBUF(1,1,1)%VOL(1:LLT)
          ENDIF
          IF (N2D == 0) THEN
            TAG22  => ELBUF_STR(NG)%GBUF%TAG22(1:LLT)
            CALL AREZO3_SIG(ALE_CONNECT,SIG,INDX,PHI,NDIM,FLUX(6*NFT+1),VOL,TAG22,IXS,ITRIMAT,NEL)
          ELSE
            NFX = NFT+(NM-1)*NUMELQ
            IF (NMULT == 0) THEN
              CALL AREZO2_SIG(ALE_CONNECT,SIG,INDX,PHI,NDIM,
     .                        FLUX(4*NFX+1),VOL,NEL)
            ELSE
              CALL BREZO2_SIG(ALE_CONNECT  ,SIG   ,IDX    ,PHI   ,NDIM  , FLUX(4*NFX+1) ,VOL,BHOLE,NM,NEL)
            ENDIF
          ENDIF
        !----------------------------!
        !        N V A R = 10        !
        !              (PLAS)        !
        !                (TB) explo. !
        !                (RK) turb.  !
        !----------------------------!
        ELSEIF (NVAR == 10 .AND. IRE == 1) THEN
           IF (MTN == 41) CYCLE
          INDX = IDX
          VOL => ELBUF_STR(NG)%BUFLY(NM)%LBUF(1,1,1)%VOL(1:LLT)
          IF (ITRIMAT /= 0 .AND. ITRIMAT /=4) THEN
            ADD0= M51_N0PHAS + (ITRIMAT-1)*M51_NVPHAS
            ADD = ADD0 + 15                                ! UVAR(ADD+1) => PLAS
            K   = LLT*(ADD-1)                               
            VAR =>ELBUF_STR(NG)%BUFLY(1)%MAT(1,1,1)%VAR(K+1:K+LLT*NDIM)
            INDX=1
          ELSEIF (MTN == 5 .OR. MTN == 97 .OR. MTN==105 .OR. ITRIMAT == 4) THEN
            VAR => ELBUF_STR(NG)%BUFLY(NM)%LBUF(1,1,1)%TB(1:LLT)
          ELSEIF (MTN == 6) THEN
            VAR => ELBUF_STR(NG)%BUFLY(NM)%LBUF(1,1,1)%RK(1:LLT)
          ELSEIF (MTN >= 28 .AND. MTN /= 67 .AND. MTN /= 49) THEN
            VAR => ELBUF_STR(NG)%BUFLY(NM)%LBUF(1,1,1)%PLA(1:LLT)
          ELSE
            VAR => ELBUF_STR(NG)%BUFLY(NM)%LBUF(1,1,1)%PLA(1:LLT)
          ENDIF        
          IF (N2D == 0) THEN 
            TAG22  => ELBUF_STR(NG)%GBUF%TAG22(1:LLT) 
            CALL AREZO3(ALE_CONNECT,VAR,IDX,PHI,NDIM,FLUX(6*NFT+1),VOL,TAG22,IXS,ITRIMAT)
          ELSE                                                        
            NFX=NFT+(NM-1)*NUMELQ                                     
            IF(NMULT == 0)THEN                                        
              CALL AREZO2(ALE_CONNECT,VAR,IDX,PHI,NDIM,FLUX(4*NFX+1),VOL)
            ELSE                                                      
              CALL BREZO2(ALE_CONNECT,VAR,IDX,PHI,NDIM,FLUX(4*NFX+1),VOL,BHOLE,NM)
            ENDIF                                                     
          ENDIF                                                       
        !----------------------------!
        !        N V A R = 11        !
        !             (EPS_P_VOL)    !
        !----------------------------!
        ELSEIF (NVAR == 11 .AND. IRE == 1 .AND. MTN == 81) THEN
          VOL => ELBUF_STR(NG)%BUFLY(NM)%LBUF(1,1,1)%VOL(1:LLT)        
          VAR => ELBUF_STR(NG)%BUFLY(NM)%MAT(1,1,1)%VAR(LLT+1:LLT*2)
          DO I=LFT,LLT                                           
            J = I+NFT                                            
            PHI(J) = VAR(I)
          ENDDO  
          IF (N2D == 0) THEN                                   
            TAG22  => ELBUF_STR(NG)%GBUF%TAG22(1:LLT)
            CALL AREZO3(ALE_CONNECT,VAR,IDX,PHI,NDIM,FLUX(6*NFT+1),VOL,TAG22,IXS,ITRIMAT)
          ELSE                                                        
            NFX=NFT+(NM-1)*NUMELQ                                     
            IF(NMULT == 0)THEN                                        
              CALL AREZO2(ALE_CONNECT,VAR,IDX,PHI,NDIM,FLUX(4*NFX+1),VOL)
            ELSE                                                      
              CALL BREZO2(ALE_CONNECT,VAR,IDX,PHI,NDIM,FLUX(4*NFX+1),VOL,BHOLE,NM)
            ENDIF                                                     
          ENDIF

        !----------------------------!
        !        N V A R = 12        !
        !              (TEMP)        !
        !----------------------------!
        ELSEIF (NVAR == 12) THEN
          INDX = IDX
          TEMP => ELBUF_STR(NG)%BUFLY(NM)%LBUF(1,1,1)%TEMP(1:LLT)
          VOL => ELBUF_STR(NG)%BUFLY(NM)%LBUF(1,1,1)%VOL(1:LLT)
          IF (ITRIMAT /= 0) THEN
            ADD0= M51_N0PHAS + (ITRIMAT-1)*M51_NVPHAS
            ADD = ADD0 + 16                                ! UVAR(ADD+1) => TEMP
            K   = LLT*(ADD-1)
            VAR =>ELBUF_STR(NG)%BUFLY(1)%MAT(1,1,1)%VAR(K+1:K+LLT*NDIM)
            INDX=1
          ELSE
            TEMP => ELBUF_STR(NG)%BUFLY(NM)%LBUF(1,1,1)%TEMP(1:LLT)
            VOL => ELBUF_STR(NG)%BUFLY(NM)%LBUF(1,1,1)%VOL(1:LLT)
          ENDIF
          IF (N2D == 0) THEN
            TAG22  => ELBUF_STR(NG)%GBUF%TAG22(1:LLT)
            CALL AREZO3(ALE_CONNECT,TEMP,IDX,PHI,NDIM,FLUX(6*NFT+1),VOL,TAG22,IXS,ITRIMAT)
          ELSE
            NFX=NFT+(NM-1)*NUMELQ
            IF(NMULT == 0)THEN
              CALL AREZO2(ALE_CONNECT,TEMP,IDX,PHI,NDIM,FLUX(4*NFX+1),VOL)
            ELSE
              CALL BREZO2(ALE_CONNECT,TEMP,IDX,PHI,NDIM,FLUX(4*NFX+1),VOL,BHOLE,NM)
            ENDIF
          ENDIF

        ENDIF !  (NVAR == 2 .AND. IRS == 1)
       ENDDO  ! next NG
      END DO !next NM
C-----------
      RETURN
      END
