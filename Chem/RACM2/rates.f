C------------------------------------------------------------------------
C
C     This routine was automatically generated by SPACK.
C
C------------------------------------------------------------------------
C------------------------------------------------------------------------
 
      subroutine rates                        (
     $    ns,nr,rk,y,w)
 
C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     This routine computes the reaction rates.
C
C------------------------------------------------------------------------
C
C     -- INPUT VARIABLES
C
C     Ns: chemical species number.
C     NR: reaction number.
C     RK: kinetic rates.
C     Y: chemical concentrations.
C
C     -- INPUT/OUTPUT VARIABLES
C
C     -- OUTPUT VARIABLES
C
C     W: reaction rates.
C
C------------------------------------------------------------------------
C
C     -- REMARKS
C
C------------------------------------------------------------------------
C
C     -- MODIFICATIONS
C
C------------------------------------------------------------------------
C
C     -- AUTHOR(S)
C
C     SPACK.
C
C------------------------------------------------------------------------
 
      implicit none
 
      integer nr,ns
      double precision rk(nr),y(ns),w(nr)
 
 
 
      w(  1) =  rk(  1) * Y(113)
      w(  2) =  rk(  2) * Y( 85)
      w(  3) =  rk(  3) * Y( 85)
      w(  4) =  rk(  4) * Y( 25)
      w(  5) =  rk(  5) * Y( 98)
      w(  6) =  rk(  6) * Y( 19)
      w(  7) =  rk(  7) * Y(116)
      w(  8) =  rk(  8) * Y(116)
      w(  9) =  rk(  9) * Y( 45)
      w( 10) =  rk( 10) * Y(117)
      w( 11) =  rk( 11) * Y(117)
      w( 12) =  rk( 12) * Y( 91)
      w( 13) =  rk( 13) * Y(106)
      w( 14) =  rk( 14) * Y( 62)
      w( 15) =  rk( 15) * Y( 31)
      w( 16) =  rk( 16) * Y(104)
      w( 17) =  rk( 17) * Y( 83)
      w( 18) =  rk( 18) * Y( 75)
      w( 19) =  rk( 19) * Y( 73)
      w( 20) =  rk( 20) * Y( 86)
      w( 21) =  rk( 21) * Y( 95)
      w( 22) =  rk( 22) * Y( 95)
      w( 23) =  rk( 23) * Y( 95)
      w( 24) =  rk( 24) * Y( 97)
      w( 25) =  rk( 25) * Y( 70)
      w( 26) =  rk( 26) * Y( 60)
      w( 27) =  rk( 27) * Y( 63)
      w( 28) =  rk( 28) * Y( 52)
      w( 29) =  rk( 29) * Y( 35)
      w( 30) =  rk( 30) * Y(112)
      w( 31) =  rk( 31) * Y( 64)
      w( 32) =  rk( 32) * Y(103)
      w( 33) =  rk( 33) * Y( 36)
      w( 34) =  rk( 34) * Y( 36)
      w( 35) =  rk( 35) * Y( 43)
      w( 36) =  rk( 36) * Y( 43) * Y( 85)
      w( 37) =  rk( 37) * Y(  1)
      w( 38) =  rk( 38) * Y(  1)
      w( 39) =  rk( 39) * Y(  1)
      w( 40) =  rk( 40) * Y( 85) * Y(114)
      w( 41) =  rk( 41) * Y( 85) * Y(115)
      w( 42) =  rk( 42) * Y(114) * Y(115)
      w( 43) =  rk( 43) * Y( 45) * Y(114)
      w( 44) =  rk( 44) * Y(115) * Y(115)
      w( 45) =  rk( 45) * Y(115) * Y(115)
      w( 46) =  rk( 46) * Y( 43) * Y(110)
      w( 47) =  rk( 47) * Y( 43) * Y(113)
      w( 48) =  rk( 48) * Y( 43) * Y(113)
      w( 49) =  rk( 49) * Y(114) * Y(110)
      w( 50) =  rk( 50) * Y(114) * Y(113)
      w( 51) =  rk( 51) * Y(114) * Y(116)
      w( 52) =  rk( 52) * Y(115) * Y(110)
      w( 53) =  rk( 53) * Y(115) * Y(113)
      w( 54) =  rk( 54) * Y( 19)
      w( 55) =  rk( 55) * Y(115) * Y(116)
      w( 56) =  rk( 56) * Y(114) * Y( 25)
      w( 57) =  rk( 57) * Y(114) * Y( 98)
      w( 58) =  rk( 58) * Y(114) * Y( 19)
      w( 59) =  rk( 59) * Y( 85) * Y(110)
      w( 60) =  rk( 60) * Y( 85) * Y(113)
      w( 61) =  rk( 61) * Y(110) * Y(110)
      w( 62) =  rk( 62) * Y(116) * Y(110)
      w( 63) =  rk( 63) * Y(116) * Y(113)
      w( 64) =  rk( 64) * Y(116) * Y(113)
      w( 65) =  rk( 65) * Y( 17)
      w( 66) =  rk( 66) * Y( 17)
      w( 67) =  rk( 67) * Y(116) * Y(116)
      w( 68) =  rk( 68) * Y(114)
      w( 69) =  rk( 69) * Y(114) * Y(  2)
      w( 70) =  rk( 70) * Y(100) * Y(114)
      w( 71) =  rk( 71) * Y( 33) * Y(114)
      w( 72) =  rk( 72) * Y( 34) * Y(114)
      w( 73) =  rk( 73) * Y( 93) * Y(114)
      w( 74) =  rk( 74) * Y( 21) * Y(114)
      w( 75) =  rk( 75) * Y( 90) * Y(114)
      w( 76) =  rk( 76) * Y(  4) * Y(114)
      w( 77) =  rk( 77) * Y(  6) * Y(114)
      w( 78) =  rk( 78) * Y(  7) * Y(114)
      w( 79) =  rk( 79) * Y(  8) * Y(114)
      w( 80) =  rk( 80) * Y( 23) * Y(114)
      w( 81) =  rk( 81) * Y( 32) * Y(114)
      w( 82) =  rk( 82) * Y( 77) * Y(114)
      w( 83) =  rk( 83) * Y( 20) * Y(114)
      w( 84) =  rk( 84) * Y( 27) * Y(114)
      w( 85) =  rk( 85) * Y( 18) * Y(114)
      w( 86) =  rk( 86) * Y( 28) * Y(114)
      w( 87) =  rk( 87) * Y( 26) * Y(114)
      w( 88) =  rk( 88) * Y( 42) * Y(114)
      w( 89) =  rk( 89) * Y( 29) * Y(114)
      w( 90) =  rk( 90) * Y( 67) * Y(114)
      w( 91) =  rk( 91) * Y( 15) * Y(114)
      w( 92) =  rk( 92) * Y( 24) * Y(114)
      w( 93) =  rk( 93) * Y(  5) * Y(114)
      w( 94) =  rk( 94) * Y(  9) * Y(114)
      w( 95) =  rk( 95) * Y( 49) * Y(110)
      w( 96) =  rk( 96) * Y( 10)
      w( 97) =  rk( 97) * Y( 50)
      w( 98) =  rk( 98) * Y( 50) * Y(110)
      w( 99) =  rk( 99) * Y( 46) * Y(110)
      w(100) =  rk(100) * Y( 13) * Y(114)
      w(101) =  rk(101) * Y( 11) * Y(114)
      w(102) =  rk(102) * Y( 51) * Y(110)
      w(103) =  rk(103) * Y( 12)
      w(104) =  rk(104) * Y( 14)
      w(105) =  rk(105) * Y( 57)
      w(106) =  rk(106) * Y( 59)
      w(107) =  rk(107) * Y( 57) * Y(110)
      w(108) =  rk(108) * Y( 58) * Y(110)
      w(109) =  rk(109) * Y( 52) * Y(114)
      w(110) =  rk(110) * Y( 53) * Y(110)
      w(111) =  rk(111) * Y( 47) * Y(110)
      w(112) =  rk(112) * Y( 48) * Y(113)
      w(113) =  rk(113) * Y( 16) * Y(114)
      w(114) =  rk(114) * Y( 55) * Y(114)
      w(115) =  rk(115) * Y( 30) * Y(114)
      w(116) =  rk(116) * Y( 30) * Y(116)
      w(117) =  rk(117) * Y( 88) * Y( 85)
      w(118) =  rk(118) * Y( 88) * Y(113)
      w(119) =  rk(119) * Y( 87) * Y(110)
      w(120) =  rk(120) * Y( 87) * Y(108)
      w(121) =  rk(121) * Y( 87) * Y(111)
      w(122) =  rk(122) * Y( 87) * Y(115)
      w(123) =  rk(123) * Y( 87) * Y(116)
      w(124) =  rk(124) * Y(117) * Y(114)
      w(125) =  rk(125) * Y( 91) * Y(114)
      w(126) =  rk(126) * Y(106) * Y(114)
      w(127) =  rk(127) * Y( 62) * Y(114)
      w(128) =  rk(128) * Y( 31) * Y(114)
      w(129) =  rk(129) * Y(104) * Y(114)
      w(130) =  rk(130) * Y( 83) * Y(114)
      w(131) =  rk(131) * Y( 95) * Y(114)
      w(132) =  rk(132) * Y( 97) * Y(114)
      w(133) =  rk(133) * Y( 75) * Y(114)
      w(134) =  rk(134) * Y( 73) * Y(114)
      w(135) =  rk(135) * Y( 86) * Y(114)
      w(136) =  rk(136) * Y( 22) * Y(114)
      w(137) =  rk(137) * Y( 60) * Y(114)
      w(138) =  rk(138) * Y( 70) * Y(114)
      w(139) =  rk(139) * Y( 63) * Y(114)
      w(140) =  rk(140) * Y( 35) * Y(114)
      w(141) =  rk(141) * Y(112) * Y(114)
      w(142) =  rk(142) * Y( 64) * Y(114)
      w(143) =  rk(143) * Y( 36) * Y(114)
      w(144) =  rk(144) * Y( 37) * Y(114)
      w(145) =  rk(145) * Y( 44) * Y(114)
      w(146) =  rk(146) * Y(103) * Y(114)
      w(147) =  rk(147) * Y( 71) * Y(114)
      w(148) =  rk(148) * Y(101) * Y(114)
      w(149) =  rk(149) * Y(117) * Y(116)
      w(150) =  rk(150) * Y( 91) * Y(116)
      w(151) =  rk(151) * Y(106) * Y(116)
      w(152) =  rk(152) * Y( 95) * Y(116)
      w(153) =  rk(153) * Y( 97) * Y(116)
      w(154) =  rk(154) * Y( 22) * Y(116)
      w(155) =  rk(155) * Y( 16) * Y(116)
      w(156) =  rk(156) * Y( 55) * Y(116)
      w(157) =  rk(157) * Y( 56) * Y(110)
      w(158) =  rk(158) * Y( 56) * Y(116)
      w(159) =  rk(159) * Y( 56) * Y(108)
      w(160) =  rk(160) * Y( 56) * Y(111)
      w(161) =  rk(161) * Y( 56) * Y(115)
      w(162) =  rk(162) * Y( 23) * Y(116)
      w(163) =  rk(163) * Y( 32) * Y(116)
      w(164) =  rk(164) * Y( 77) * Y(116)
      w(165) =  rk(165) * Y( 20) * Y(116)
      w(166) =  rk(166) * Y( 27) * Y(116)
      w(167) =  rk(167) * Y( 75) * Y(116)
      w(168) =  rk(168) * Y( 86) * Y(116)
      w(169) =  rk(169) * Y( 18) * Y(116)
      w(170) =  rk(170) * Y( 28) * Y(116)
      w(171) =  rk(171) * Y( 44) * Y(116)
      w(172) =  rk(172) * Y( 23) * Y( 85)
      w(173) =  rk(173) * Y( 32) * Y( 85)
      w(174) =  rk(174) * Y( 77) * Y( 85)
      w(175) =  rk(175) * Y( 20) * Y( 85)
      w(176) =  rk(176) * Y( 27) * Y( 85)
      w(177) =  rk(177) * Y( 75) * Y( 85)
      w(178) =  rk(178) * Y( 73) * Y( 85)
      w(179) =  rk(179) * Y( 86) * Y( 85)
      w(180) =  rk(180) * Y( 18) * Y( 85)
      w(181) =  rk(181) * Y( 28) * Y( 85)
      w(182) =  rk(182) * Y( 22) * Y( 85)
      w(183) =  rk(183) * Y( 63) * Y( 85)
      w(184) =  rk(184) * Y( 54) * Y(113)
      w(185) =  rk(185) * Y(111) * Y(113)
      w(186) =  rk(186) * Y(102) * Y(113)
      w(187) =  rk(187) * Y( 36)
      w(188) =  rk(188) * Y( 37)
      w(189) =  rk(189) * Y( 99) * Y(113)
      w(190) =  rk(190) * Y( 44)
      w(191) =  rk(191) * Y(108) * Y(110)
      w(192) =  rk(192) * Y(107) * Y(110)
      w(193) =  rk(193) * Y(105) * Y(110)
      w(194) =  rk(194) * Y( 61) * Y(110)
      w(195) =  rk(195) * Y( 81) * Y(110)
      w(196) =  rk(196) * Y( 84) * Y(110)
      w(197) =  rk(197) * Y( 69) * Y(110)
      w(198) =  rk(198) * Y( 78) * Y(110)
      w(199) =  rk(199) * Y( 72) * Y(110)
      w(200) =  rk(200) * Y( 99) * Y(110)
      w(201) =  rk(201) * Y( 74) * Y(110)
      w(202) =  rk(202) * Y( 68) * Y(110)
      w(203) =  rk(203) * Y( 89) * Y(110)
      w(204) =  rk(204) * Y( 66) * Y(110)
      w(205) =  rk(205) * Y( 65) * Y(110)
      w(206) =  rk(206) * Y( 76) * Y(110)
      w(207) =  rk(207) * Y( 59) * Y(110)
      w(208) =  rk(208) * Y( 82) * Y(110)
      w(209) =  rk(209) * Y(111) * Y(110)
      w(210) =  rk(210) * Y(102) * Y(110)
      w(211) =  rk(211) * Y( 96) * Y(110)
      w(212) =  rk(212) * Y( 92) * Y(110)
      w(213) =  rk(213) * Y( 80) * Y(110)
      w(214) =  rk(214) * Y( 79) * Y(110)
      w(215) =  rk(215) * Y(108) * Y(115)
      w(216) =  rk(216) * Y(107) * Y(115)
      w(217) =  rk(217) * Y(105) * Y(115)
      w(218) =  rk(218) * Y( 61) * Y(115)
      w(219) =  rk(219) * Y( 81) * Y(115)
      w(220) =  rk(220) * Y( 84) * Y(115)
      w(221) =  rk(221) * Y( 69) * Y(115)
      w(222) =  rk(222) * Y( 78) * Y(115)
      w(223) =  rk(223) * Y( 72) * Y(115)
      w(224) =  rk(224) * Y( 99) * Y(115)
      w(225) =  rk(225) * Y( 74) * Y(115)
      w(226) =  rk(226) * Y( 89) * Y(115)
      w(227) =  rk(227) * Y( 66) * Y(115)
      w(228) =  rk(228) * Y( 65) * Y(115)
      w(229) =  rk(229) * Y( 76) * Y(115)
      w(230) =  rk(230) * Y( 50) * Y(115)
      w(231) =  rk(231) * Y( 57) * Y(115)
      w(232) =  rk(232) * Y( 59) * Y(115)
      w(233) =  rk(233) * Y( 82) * Y(115)
      w(234) =  rk(234) * Y( 54) * Y(115)
      w(235) =  rk(235) * Y(111) * Y(115)
      w(236) =  rk(236) * Y(102) * Y(115)
      w(237) =  rk(237) * Y( 92) * Y(115)
      w(238) =  rk(238) * Y( 68) * Y(115)
      w(239) =  rk(239) * Y( 96) * Y(115)
      w(240) =  rk(240) * Y( 80) * Y(115)
      w(241) =  rk(241) * Y( 79) * Y(115)
      w(242) =  rk(242) * Y(108) * Y(108)
      w(243) =  rk(243) * Y(107) * Y(108)
      w(244) =  rk(244) * Y(105) * Y(108)
      w(245) =  rk(245) * Y( 61) * Y(108)
      w(246) =  rk(246) * Y( 81) * Y(108)
      w(247) =  rk(247) * Y( 84) * Y(108)
      w(248) =  rk(248) * Y( 69) * Y(108)
      w(249) =  rk(249) * Y( 78) * Y(108)
      w(250) =  rk(250) * Y( 72) * Y(108)
      w(251) =  rk(251) * Y( 99) * Y(108)
      w(252) =  rk(252) * Y( 74) * Y(108)
      w(253) =  rk(253) * Y( 89) * Y(108)
      w(254) =  rk(254) * Y( 66) * Y(108)
      w(255) =  rk(255) * Y( 65) * Y(108)
      w(256) =  rk(256) * Y( 76) * Y(108)
      w(257) =  rk(257) * Y( 50) * Y(108)
      w(258) =  rk(258) * Y( 49) * Y(108)
      w(259) =  rk(259) * Y( 46) * Y(108)
      w(260) =  rk(260) * Y( 59) * Y(108)
      w(261) =  rk(261) * Y( 51) * Y(108)
      w(262) =  rk(262) * Y( 57) * Y(108)
      w(263) =  rk(263) * Y( 58) * Y(108)
      w(264) =  rk(264) * Y( 53) * Y(108)
      w(265) =  rk(265) * Y( 47) * Y(108)
      w(266) =  rk(266) * Y( 82) * Y(108)
      w(267) =  rk(267) * Y(111) * Y(108)
      w(268) =  rk(268) * Y(102) * Y(108)
      w(269) =  rk(269) * Y( 92) * Y(108)
      w(270) =  rk(270) * Y( 68) * Y(108)
      w(271) =  rk(271) * Y( 96) * Y(108)
      w(272) =  rk(272) * Y( 80) * Y(108)
      w(273) =  rk(273) * Y( 79) * Y(108)
      w(274) =  rk(274) * Y(107) * Y(111)
      w(275) =  rk(275) * Y(105) * Y(111)
      w(276) =  rk(276) * Y( 61) * Y(111)
      w(277) =  rk(277) * Y( 81) * Y(111)
      w(278) =  rk(278) * Y( 84) * Y(111)
      w(279) =  rk(279) * Y( 69) * Y(111)
      w(280) =  rk(280) * Y( 78) * Y(111)
      w(281) =  rk(281) * Y( 72) * Y(111)
      w(282) =  rk(282) * Y( 99) * Y(111)
      w(283) =  rk(283) * Y( 74) * Y(111)
      w(284) =  rk(284) * Y( 89) * Y(111)
      w(285) =  rk(285) * Y( 66) * Y(111)
      w(286) =  rk(286) * Y( 65) * Y(111)
      w(287) =  rk(287) * Y( 76) * Y(111)
      w(288) =  rk(288) * Y( 50) * Y(111)
      w(289) =  rk(289) * Y( 49) * Y(111)
      w(290) =  rk(290) * Y( 46) * Y(111)
      w(291) =  rk(291) * Y( 58) * Y(111)
      w(292) =  rk(292) * Y( 53) * Y(111)
      w(293) =  rk(293) * Y( 47) * Y(111)
      w(294) =  rk(294) * Y( 59) * Y(111)
      w(295) =  rk(295) * Y( 51) * Y(111)
      w(296) =  rk(296) * Y( 57) * Y(111)
      w(297) =  rk(297) * Y( 82) * Y(111)
      w(298) =  rk(298) * Y(111) * Y(111)
      w(299) =  rk(299) * Y(102) * Y(102)
      w(300) =  rk(300) * Y( 96) * Y(111)
      w(301) =  rk(301) * Y( 92) * Y(111)
      w(302) =  rk(302) * Y( 68) * Y(111)
      w(303) =  rk(303) * Y( 80) * Y(111)
      w(304) =  rk(304) * Y( 79) * Y(111)
      w(305) =  rk(305) * Y( 80) * Y( 80)
      w(306) =  rk(306) * Y( 80) * Y( 79)
      w(307) =  rk(307) * Y( 79) * Y( 79)
      w(308) =  rk(308) * Y(108) * Y(116)
      w(309) =  rk(309) * Y(107) * Y(116)
      w(310) =  rk(310) * Y(105) * Y(116)
      w(311) =  rk(311) * Y( 61) * Y(116)
      w(312) =  rk(312) * Y( 81) * Y(116)
      w(313) =  rk(313) * Y( 84) * Y(116)
      w(314) =  rk(314) * Y( 69) * Y(116)
      w(315) =  rk(315) * Y( 78) * Y(116)
      w(316) =  rk(316) * Y( 72) * Y(116)
      w(317) =  rk(317) * Y( 66) * Y(116)
      w(318) =  rk(318) * Y( 65) * Y(116)
      w(319) =  rk(319) * Y( 76) * Y(116)
      w(320) =  rk(320) * Y( 59) * Y(116)
      w(321) =  rk(321) * Y( 82) * Y(116)
      w(322) =  rk(322) * Y(111) * Y(116)
      w(323) =  rk(323) * Y(102) * Y(116)
      w(324) =  rk(324) * Y( 68) * Y(116)
      w(325) =  rk(325) * Y( 96) * Y(116)
      w(326) =  rk(326) * Y( 92) * Y(116)
      w(327) =  rk(327) * Y( 80) * Y(116)
      w(328) =  rk(328) * Y( 79) * Y(116)
      w(329) =  rk(329) * Y( 74) * Y(116)
      w(330) =  rk(330) * Y( 99) * Y(116)
      w(331) =  rk(331) * Y( 49) * Y(116)
      w(332) =  rk(332) * Y( 50) * Y(116)
      w(333) =  rk(333) * Y( 46) * Y(116)
      w(334) =  rk(334) * Y( 51) * Y(116)
      w(335) =  rk(335) * Y( 57) * Y(116)
      w(336) =  rk(336) * Y( 58) * Y(116)
      w(337) =  rk(337) * Y( 53) * Y(116)
      w(338) =  rk(338) * Y( 47) * Y(116)
      w(339) =  rk(339) * Y( 89) * Y(116)
      w(340) =  rk(340) * Y(109) * Y(115)
      w(341) =  rk(341) * Y(109) * Y(108)
      w(342) =  rk(342) * Y(109) * Y(111)
      w(343) =  rk(343) * Y(109) * Y(109)
      w(344) =  rk(344) * Y(109) * Y(110)
      w(345) =  rk(345) * Y(109) * Y(116)
      w(346) =  rk(346) * Y( 94) * Y(110)
      w(347) =  rk(347) * Y( 94) * Y(115)
      w(348) =  rk(348) * Y( 94) * Y(108)
      w(349) =  rk(349) * Y( 94) * Y(111)
      w(350) =  rk(350) * Y( 94) * Y(116)
      w(351) =  rk(351) * Y(102) * Y(111)
      w(352) =  rk(352) * Y( 38)
      w(353) =  rk(353) * Y( 38) * Y( 40)
      w(354) =  rk(354) * Y( 39) * Y( 40)
      w(355) =  rk(355) * Y( 32)
      w(356) =  rk(356) * Y( 77)
      w(357) =  rk(357) * Y(  9)
      w(358) =  rk(358) * Y( 11)
      w(359) =  rk(359) * Y(117)
      w(360) =  rk(360) * Y(106)
      w(361) =  rk(361) * Y(104)
      w(362) =  rk(362) * Y( 27)
      w(363) =  rk(363) * Y( 55)
      w(364) =  rk(364) * Y( 71)
      w(365) =  rk(365) * Y(101)
      w(366) =  rk(366) * Y( 75)
      w(367) =  rk(367) * Y( 93)
      w(368) =  rk(368) * Y( 91)
      w(369) =  rk(369) * Y( 62)
      w(370) =  rk(370) * Y(  5)
      w(371) =  rk(371) * Y( 16)
      w(372) =  rk(372) * Y( 52)
      w(373) =  rk(373) * Y( 90)
      w(374) =  rk(374) * Y( 86)
      w(375) =  rk(375) * Y( 28)
      w(376) =  rk(376) * Y(110)
      w(377) =  rk(377) * Y(113)
      w(378) =  rk(378) * Y( 41)
      w(379) =  rk(379) * Y( 25)
      w(380) =  rk(380) * Y( 85)
      w(381) =  rk(381) * Y(114)
      w(382) =  rk(382) * Y(115)
      w(383) =  rk(383) * Y( 45)
      w(384) =  rk(384) * Y( 98)
      w(385) =  rk(385) * Y( 19)
      w(386) =  rk(386) * Y(116)
      w(387) =  rk(387) * Y( 38)
      w(388) =  rk(388) * Y( 39)
      w(389) =  rk(389) * Y( 41)
 
      RETURN
      END
 
