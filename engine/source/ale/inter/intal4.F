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
      !||    intal4          ../engine/source/ale/inter/intal4.F
      !||--- called by ------------------------------------------------------
      !||    intal3          ../engine/source/ale/inter/intal3.F
      !||--- calls      -----------------------------------------------------
      !||    iqela3          ../engine/source/ale/inter/iqela3.F
      !||--- uses       -----------------------------------------------------
      !||====================================================================
      SUBROUTINE INTAL4(IPARI,INTBUF_TAB         ,X         ,V         ,
     1           A           ,VG      ,ICODE     ,IBCS      ,BCS       )
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
       USE INTBUFDEF_MOD
C-----------------------------------------------
C   I m p l i c i t   T y p e s
C-----------------------------------------------
#include      "implicit_f.inc"
C-----------------------------------------------
C   C o m m o n   B l o c k s
C-----------------------------------------------
#include      "com04_c.inc"
#include      "tabsiz_c.inc"
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      my_real X(3,NUMNOD), V(3,NUMNOD), A(3,NUMNOD), VG(*), BCS(*)      
      INTEGER IPARI(SIPARI), ICODE(*), IBCS(*)
      TYPE(INTBUF_STRUCT_) INTBUF_TAB
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: NSN 
C-----------------------------------------------
C   S o u r c e   L i n e s
C-----------------------------------------------
      NSN=IPARI(5)
      CALL IQELA3
     1  (V               ,VG              ,INTBUF_TAB%IRECTM,INTBUF_TAB%CSTS,INTBUF_TAB%NSV,
     2   INTBUF_TAB%ILOCS,INTBUF_TAB%IRTLM,ICODE            ,IBCS           ,BCS       ,
     3   INTBUF_TAB%MSR, NSN)
C
      RETURN
      END
