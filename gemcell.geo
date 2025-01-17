//
//
//
//
thDielec  = 0.0680;  // dielectric thickness
//thDielec  = #####THICK#####;  // dielectric thickness
thCopper  = 0.00051;  // copper thickness
radCenter = 0.0085;  // the hole radius
radSurfac = 0.0085;  // the hole radius 
//radCenter = #####PHI#####;  // the hole radius
//radSurfac = #####PHI#####;  // the hole radius 
pitchHoles= 0.0280;  // the "pitch", or distance between GEM holes

distUpperElec = 0.800;  // distance from GEM plates to upper exterior electrode
distLowerPads = 0.800;  // distance from lower LEM plate to pad (readout) plane

//lcDielectricHole   = 0.008;
//lcEtchingHole      = 0.006;
//lcCopperPlateBdry  = 0.006;
//lcExtElectrodeBdry = 0.0400;
//lcGEMHole          = 0.015;

lcDielectricHole   = 0.0008;
lcEtchingHole      = 0.0006;
lcCopperPlateBdry  = 0.0006;
lcExtElectrodeBdry = 0.0400;
lcGEMHole          = 0.0015;

// *******************************************************
// Hole 1 (quarter hole)
// *******************************************************

// ----- Center
Point(0) = {0, 0,  0,         lcGEMHole};
Point(1) = {0, 0,  thDielec/2,lcGEMHole};
Point(2) = {0, 0, -thDielec/2,lcGEMHole};
Point(3) = {0, 0,  (thDielec/2+thCopper),lcGEMHole};
Point(4) = {0, 0, -(thDielec/2+thCopper),lcGEMHole};

// ----- Dielectric hole
// top
Point(5) = {radSurfac, 0, thDielec/2, lcDielectricHole};
Point(6) = {0, radSurfac, thDielec/2, lcDielectricHole};
Circle(200) = {5, 1, 6};
// center
Point(7) = {radCenter, 0, 0, lcDielectricHole};
Point(8) = {0, radCenter, 0, lcDielectricHole};
Circle(201) = {7, 0, 8};
// bottom
Point(9) = {radSurfac, 0, -thDielec/2, lcDielectricHole};
Point(10)= {0, radSurfac, -thDielec/2, lcDielectricHole};
Circle(202) = {9, 2, 10};



// ----- Upper Etching
// top
Point(11) = {radSurfac, 0, thDielec/2+thCopper,lcEtchingHole};
Point(12) = {0, radSurfac, thDielec/2+thCopper,lcEtchingHole};
// bottom
//Point(13) = {radSurfac, 0, thDielec/2,lcEtchingHole};
//Point(14) = {0, radSurfac, thDielec/2,lcEtchingHole};

// circular boundary
Circle(203) = {11, 3, 12};
//Circle(204) = {13, 1, 14};

// ----- Lower Etching
// top
//Point(15) = {radSurfac, 0, -thDielec/2,lcEtchingHole};
//Point(16) = {0, radSurfac, -thDielec/2,lcEtchingHole};

// bottom
Point(17) = {radSurfac, 0, -(thDielec/2+thCopper),lcEtchingHole};
Point(18) = {0, radSurfac, -(thDielec/2+thCopper),lcEtchingHole};

// circular boundaries
//Circle(205) = {15, 2, 16};
Circle(206) = {17, 4, 18};

// Lines connecting top and bottom
Line(207) = {9, 7};
Line(208) = {7, 5};
Line(209) = {10,8};
Line(210) = {8, 6};


// *******************************************************
// Hole 2 (quarter hole)
// *******************************************************

// ----- Center
Point(20) = {pitchHoles/2, pitchHoles*Sqrt(3)/2,  0,         lcGEMHole};
Point(21) = {pitchHoles/2, pitchHoles*Sqrt(3)/2,  thDielec/2,lcGEMHole};
Point(22) = {pitchHoles/2, pitchHoles*Sqrt(3)/2, -thDielec/2,lcGEMHole};
Point(23) = {pitchHoles/2, pitchHoles*Sqrt(3)/2,  (thDielec/2+thCopper),lcGEMHole};
Point(24) = {pitchHoles/2, pitchHoles*Sqrt(3)/2, -(thDielec/2+thCopper),lcGEMHole};

// ----- Dielectric hole
// top
Point(25) = {pitchHoles/2-radSurfac, pitchHoles*Sqrt(3)/2,           thDielec/2, lcDielectricHole};
Point(26) = {pitchHoles/2,           pitchHoles*Sqrt(3)/2-radSurfac, thDielec/2, lcDielectricHole};
Circle(220) = {25, 21, 26};
// center
Point(27) = {pitchHoles/2-radCenter, pitchHoles*Sqrt(3)/2,           0, lcDielectricHole};
Point(28) = {pitchHoles/2,           pitchHoles*Sqrt(3)/2-radCenter, 0, lcDielectricHole};
Circle(221) = {27, 20, 28};
// bottom
Point(29) = {pitchHoles/2-radSurfac, pitchHoles*Sqrt(3)/2,           -thDielec/2, lcDielectricHole};
Point(30) = {pitchHoles/2,           pitchHoles*Sqrt(3)/2-radSurfac, -thDielec/2, lcDielectricHole};
Circle(222) = {29, 22, 30};


// ----- Upper Etching
// top
Point(31) = {pitchHoles/2-radSurfac, pitchHoles*Sqrt(3)/2,           thDielec/2+thCopper,lcEtchingHole};
Point(32) = {pitchHoles/2,           pitchHoles*Sqrt(3)/2-radSurfac, thDielec/2+thCopper,lcEtchingHole};
// bottom
//Point(33) = {pitchHoles/2-radSurfac, pitchHoles*Sqrt(3)/2,           thDielec/2,lcEtchingHole};
//Point(34) = {pitchHoles/2,           pitchHoles*Sqrt(3)/2-radSurfac, thDielec/2,lcEtchingHole};

// circular boundary
Circle(223) = {31, 23, 32};
//Circle(224) = {33, 21, 34};

// ----- Lower Etching
// top
//Point(35) = {pitchHoles/2-radSurfac, pitchHoles*Sqrt(3)/2,           -thDielec/2,lcEtchingHole};
//Point(36) = {pitchHoles/2,           pitchHoles*Sqrt(3)/2-radSurfac, -thDielec/2,lcEtchingHole};

// bottom
Point(37) = {pitchHoles/2-radSurfac, pitchHoles*Sqrt(3)/2,          -(thDielec/2+thCopper),lcEtchingHole};
Point(38) = {pitchHoles/2,           pitchHoles*Sqrt(3)/2-radSurfac,-(thDielec/2+thCopper),lcEtchingHole};

// circular boundaries
//Circle(225) = {35, 22, 36};
Circle(226) = {37, 24, 38};

// Lines connecting top and bottom
Line(227) = {29, 27};
Line(228) = {27, 25};
Line(229) = {30, 28};
Line(230) = {28, 26};


// *******************************************************
// Copper planes
// *******************************************************

// ----- Upper copper
// points for corners
Point(40) = {pitchHoles/2, 0, thDielec/2,          lcCopperPlateBdry};
Point(41) = {pitchHoles/2, 0, thDielec/2+thCopper, lcCopperPlateBdry};

Point(42) = {0, pitchHoles*Sqrt(3)/2, thDielec/2,          lcCopperPlateBdry};
Point(43) = {0, pitchHoles*Sqrt(3)/2, thDielec/2+thCopper, lcCopperPlateBdry};

// make upper boundary
Line(240) = {11,  41};
Line(241) = {41, 32};
Line(242) = {31,43};
Line(243) = {43,12};
// make lower boundary
Line(244) = {5,  40};
Line(245) = {40, 26};
Line(246) = {25,42};
Line(247) = {42,6};
// Connect the upper and lower points with lines to form the plate
Line(248) = {5,11};
Line(249) = {40,41};
Line(250) = {26,32};
Line(251) = {25,31};
Line(252) = {42,43};
Line(253) = {6,12};

// ----- Lower copper
// points for corners
Point(44) = {pitchHoles/2, 0, -thDielec/2,          lcCopperPlateBdry};
Point(45) = {pitchHoles/2, 0, -(thDielec/2+thCopper), lcCopperPlateBdry};

Point(46) = {0, pitchHoles*Sqrt(3)/2, -thDielec/2,          lcCopperPlateBdry};
Point(47) = {0, pitchHoles*Sqrt(3)/2, -(thDielec/2+thCopper), lcCopperPlateBdry};

// make upper boundary
Line(254) = {9,  44};
Line(255) = {44, 30};
Line(256) = {29,46};
Line(257) = {46,10};
// make lower boundary
Line(258) = {17,  45};
Line(259) = {45, 38};
Line(260) = {37,47};
Line(261) = {47,18};
// Connect the upper and lower points with lines to form the plate
Line(262) = {17,9};
Line(263) = {45,44};
Line(264) = {38,30};
Line(265) = {37,29};
Line(266) = {47,46};
Line(267) = {18,10};


// **********************************************
// External Electrodes
// **********************************************

// Top electrode
Point(50) = {0,            0, thDielec/2+thCopper+distUpperElec,lcExtElectrodeBdry};
Point(51) = {pitchHoles/2, 0, thDielec/2+thCopper+distUpperElec,lcExtElectrodeBdry};
Point(52) = {pitchHoles/2, pitchHoles*Sqrt(3)/2, thDielec/2+thCopper+distUpperElec,lcExtElectrodeBdry};
Point(53) = {0,            pitchHoles*Sqrt(3)/2, thDielec/2+thCopper+distUpperElec,lcExtElectrodeBdry};

// Top electrode lines
Line(270) = {50,51};
Line(271) = {51,52};
Line(272) = {52,53};
Line(273) = {53,50};

// Connect the top electrode to the LEM.
Line(274) = {41,51};
Line(275) = { 3,50};
Line(276) = {23,52};
Line(277) = {43,53};

// Bottom electrode
Point(60) = {0,            0, -(thDielec/2+thCopper)-distLowerPads,lcExtElectrodeBdry};
Point(61) = {pitchHoles/2, 0, -(thDielec/2+thCopper)-distLowerPads,lcExtElectrodeBdry};
Point(62) = {pitchHoles/2, pitchHoles*Sqrt(3)/2, -(thDielec/2+thCopper)-distLowerPads,lcExtElectrodeBdry};
Point(63) = {0,            pitchHoles*Sqrt(3)/2, -(thDielec/2+thCopper)-distLowerPads,lcExtElectrodeBdry};

// Bottom electrode lines
Line(280) = {60,61};
Line(281) = {61,62};
Line(282) = {62,63};
Line(283) = {63,60};

// Connect the bottom electrode to the LEM.
Line(284) = {61,45};
Line(285) = {60,4};
Line(286) = {62,24};
Line(287) = {63,47};

Line(288) = {40,44};
Line(289) = {46,42};

Line(300) = {3,11};
Line(301) = {3,12};
Line(302) = {23,31};
Line(303) = {23,32};
Line(304) = {4,18};
Line(305) = {4,17};
Line(306) = {24,38};
Line(307) = {24,37};

Line(308) = {1,3};
Line(309) = {1,6};
Line(310) = {1,5};
Line(311) = {21,23};
Line(312) = {21,26};
Line(313) = {21,25};
Line(314) = {2,4};
Line(315) = {2,9};
Line(316) = {2,10};
Line(317) = {22,24};
Line(318) = {22,29};
Line(319) = {22,30};

Line(320) = {0,1};
Line(321) = {0,8};
Line(322) = {0,7};
Line(323) = {0,2};
Line(324) = {20,21};
Line(325) = {20,22};
Line(326) = {20,27};
Line(327) = {20,28};


// *************************************************
// Define surfaces
// *************************************************

// Copper plate surfacese upper
Line Loop(400) = {253, -203, -248, 200};
Line Loop(401) = {250, -223, -251, 220};

Ruled Surface(402) = {400};
Ruled Surface(403) = {401};

Line Loop(404) = {240, -249, -244,  248};
Line Loop(405) = {249,  241, -250, -245};
Line Loop(406) = {251,  242, -252, -246};
Line Loop(407) = {252,  243, -253, -247};
Line Loop(408) = {203, -243, -242, 223, -241, -240};
Line Loop(409) = {200, -247, -246, 220, -245, -244};

Plane Surface(410) = {404};
Plane Surface(411) = {405};
Plane Surface(412) = {406};
Plane Surface(413) = {407};
Plane Surface(414) = {408};
Plane Surface(415) = {409};

// Dielectric surfaces
Line Loop(416) = {201,  210, -200, -208};
Line Loop(417) = {201, -209, -202, 207};
Line Loop(418) = {221,  230, -220, -228};
Line Loop(419) = {221, -229, -222, 227};

Ruled Surface(420) = {416};
Ruled Surface(421) = {417};
Ruled Surface(422) = {418};
Ruled Surface(423) = {419};

Line Loop(424) = {207, 208, 244,  288,-254};
Line Loop(425) = {-288, 245, -230, -229, -255};
Line Loop(426) = {228,  246, -289, -256, 227};
Line Loop(427) = {289,  247, -210, -209, -257};

//Line Loop(428) = {203, -243, -242, 223, -241, -240};
//Line Loop(429) = {200, -247, -246, 220, -245, -244};

Plane Surface(430) = {424};
Plane Surface(431) = {425};
Plane Surface(432) = {426};
Plane Surface(433) = {427};
//Plane Surface(434) = {428}; // same with 414
//Plane Surface(435) = {429}; // same with 415


// Copper plate surfacese lower
Line Loop(450) = {-262, 206, 267, -202};
Line Loop(451) = {264, -222, -265, 226};

Ruled Surface(452) = {450};
Ruled Surface(453) = {451};

Line Loop(454) = {-258, 262, 254,  -263};
Line Loop(455) = {263,  255, -264, -259};
Line Loop(456) = {-260,  265, 256, -266};
Line Loop(457) = {-261,  266, 257, -267};
Line Loop(458) = {254, 255, -222, 256, 257, -202};
Line Loop(459) = {258, 259, -226, 260, 261, -206};

Plane Surface(460) = {454};
Plane Surface(461) = {455};
Plane Surface(462) = {456};
Plane Surface(463) = {457};
Plane Surface(464) = {458};
Plane Surface(465) = {459};





// Bounding surfaces

// lower electoro
Line Loop(470) = {280,281,282,283};
// upper electro
Line Loop(471) = {270,271, 272, 273};

Plane Surface(472) = {470};
Plane Surface(473) = {471};

// upper drift volum
Line Loop(474) = {275,270,-274,-240, -300};
Line Loop(475) = {274,271,-276,303, -241};
Line Loop(476) = {276,272,-277, -242, -302};
Line Loop(477) = {277, 273, -275,301,-243};

Line Loop(478) = {311,302,-251,-313};
Line Loop(479) = {250, -303, -311, 312};
Line Loop(480) = {253, -301, -308, 309};
Line Loop(481) = {308, 300,-248, -310};

Line Loop(482) = {320, 310, -208, -322};
Line Loop(483) = {210, -309, -320, 321};
Line Loop(484) = {230, -312, -324, 327};
Line Loop(485) = {324, 313,-228, -326};

Plane Surface(486) = {474};
//Plane Surface(487) = {475};
Plane Surface(488) = {475};
Plane Surface(489) = {476};
Plane Surface(490) = {477};
Plane Surface(491) = {478};
Plane Surface(492) = {479};
Plane Surface(493) = {480};
Plane Surface(494) = {481};
Plane Surface(495) = {482};
Plane Surface(496) = {483};
Plane Surface(497) = {484};
Plane Surface(498) = {485};

// lower drift volum
Line Loop(500) = {-323,322,-207, -315};
Line Loop(501) = {209, -321,323, 316};
Line Loop(502) = {229, -327, 325, 319};
Line Loop(503) = {-325, 326,-227, -318};

Line Loop(504) = {264, -319, 317, 306};
Line Loop(505) = {-317, 318, -265, -307};
Line Loop(506) = {267, -316, 314, 304};
Line Loop(507) = {-314, 315, -262, -305};

Line Loop(508) = {-280, 285, 305, 258, -284};
Line Loop(509) = {284,259,-306,-286,-281 };
Line Loop(510) = {286, 307,260, -287, -282};
Line Loop(511) = {287, 261,-304, -285, -283};

Plane Surface(512) = {500};
Plane Surface(513) = {501};
Plane Surface(514) = {502};
Plane Surface(515) = {503};
Plane Surface(516) = {504};
Plane Surface(517) = {505};
Plane Surface(518) = {506};
Plane Surface(519) = {507};
Plane Surface(520) = {508};
Plane Surface(521) = {509};
Plane Surface(522) = {510};
Plane Surface(523) = {511};


// Volumes
// dielectric
Surface Loop(600) = {420,421,430,431,422,423,432,433,415,464};
Volume(601) = {600};
// upper coppe
Surface Loop(602) = {402,410,411,403,412,413,415,414};
Volume(603) = {602};
// lower coppe
Surface Loop(604) = {452,460,461,453,462,463,464,465};
Volume(605) = {604};


// gas
//Surface Loop(606) = {472,473,  489,491,412,432,498,515,517,462,522,  490,413,493,433,496,513,463,518,523,  486,494,410,495,430,512,519,460,520,  488,411,492,431,497,514,461,516,521};
//Surface Loop(606) = {472,473,  489,491,412,432,498,515,517,462,522,  490,413,493,433,496,513,463,518,523,  486,494,410,495,430,512,519,460,520,  488,411,492,431,497,514,461,516,521,  -413,-433,-463,-402,-420,-421,-452,-410,-430,-460,-411,-431,-461,-403,-422,-423,-453,-412,-432,-462,-414,-465};
Surface Loop(606) = {472,473,  489,491,498,515,517,522,  490,493,496,513,518,523,  486,494,495,512,519,520,  488,492,497,514,516,521,  -402,-420,-421,-452,-403,-422,-423,-453,-414,-465};
Volume(607) = {606};


// Physical surfaces

// Surfaces to which voltages will be applied
//Physical Surface(physsurf_upper_cp) = {
Physical Surface(700) = {402,410,411,403,412,413,415,414};

//Physical Surface(physsurf_lower_cp) = {
Physical Surface(701) = {452,460,461,453,462,463,464,465};

//Physical Surface(physsurf_upper_el) = {ps_bsurf9};
//Physical Surface(physsurf_lower_el) = {ps_bsurf14};
Physical Surface(702) = {473};
Physical Surface(703) = {472};

// Surfaces for periodic boundary conditions
//Physical Surface(704) = {486, 494,410, 430, 495,512, 519, 460, 520};
//Physical Surface(704) = {486, 494, 430, 495,512, 519, 520};
//Physical Surface(705) = {488, 411, 492, 431, 497, 514, 516, 461, 521};
//Physical Surface(705) = {488, 492, 431, 497, 514, 516, 521};

//Physical Surface(706) = {489,491, 412, 498, 432, 515, 517, 462, 522};
//Physical Surface(706) = {489,491, 498, 432, 515, 517, 522};
//Physical Surface(707) = {490, 413, 493, 496, 433, 513, 518, 523, 463};
//Physical Surface(707) = {490, 493, 496, 433, 513, 518, 523};


// Physical volumes
//Physical Volume(physvol_gas) = {vol_gas};
//Physical Volume(physvol_dielectric) = {vol_dielectric};
//Physical Volume(physvol_upper_cp) = {vol_upper_cp};
//Physical Volume(physvol_lower_cp) = {vol_lower_cp};

Physical Volume(704) = {607};
Physical Volume(705) = {601};
Physical Volume(706) = {603};
Physical Volume(707) = {605};













