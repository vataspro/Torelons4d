module operator_dispatch
    use iso_fortran_env, only : real64
    use parameters
    use lattice
    use thermal_lines
    implicit none

    integer, parameter :: NUMBER_OF_LOOPS = 337
    integer, parameter :: LOOP_NUMBERS(NUMBER_OF_LOOPS) = [ &
        1, &
        2, &
        3, &
        4, &
        5, &
        6, &
        7, &
        8, &
        9, &
        10, &
        11, &
        12, &
        13, &
        14, &
        15, &
        16, &
        17, &
        18, &
        19, &
        20, &
        21, &
        22, &
        23, &
        24, &
        25, &
        26, &
        27, &
        28, &
        29, &
        30, &
        31, &
        32, &
        33, &
        34, &
        35, &
        36, &
        37, &
        38, &
        39, &
        40, &
        41, &
        42, &
        43, &
        44, &
        45, &
        46, &
        47, &
        48, &
        49, &
        50, &
        51, &
        52, &
        53, &
        54, &
        55, &
        56, &
        57, &
        58, &
        59, &
        60, &
        61, &
        62, &
        63, &
        64, &
        65, &
        66, &
        67, &
        68, &
        69, &
        70, &
        71, &
        72, &
        73, &
        74, &
        75, &
        76, &
        77, &
        78, &
        79, &
        80, &
        81, &
        82, &
        83, &
        84, &
        85, &
        86, &
        87, &
        88, &
        89, &
        90, &
        91, &
        92, &
        93, &
        94, &
        95, &
        96, &
        97, &
        98, &
        99, &
        100, &
        101, &
        102, &
        103, &
        104, &
        105, &
        106, &
        107, &
        108, &
        109, &
        110, &
        111, &
        112, &
        113, &
        114, &
        115, &
        116, &
        117, &
        118, &
        119, &
        120, &
        121, &
        122, &
        123, &
        124, &
        125, &
        126, &
        127, &
        128, &
        129, &
        130, &
        131, &
        132, &
        133, &
        134, &
        135, &
        136, &
        137, &
        138, &
        139, &
        140, &
        141, &
        142, &
        143, &
        144, &
        145, &
        146, &
        147, &
        148, &
        149, &
        150, &
        151, &
        152, &
        153, &
        154, &
        155, &
        156, &
        157, &
        158, &
        159, &
        160, &
        161, &
        162, &
        163, &
        164, &
        165, &
        166, &
        167, &
        168, &
        169, &
        170, &
        171, &
        172, &
        173, &
        174, &
        175, &
        176, &
        177, &
        178, &
        179, &
        180, &
        181, &
        182, &
        183, &
        184, &
        185, &
        186, &
        187, &
        188, &
        189, &
        190, &
        191, &
        192, &
        193, &
        194, &
        195, &
        196, &
        197, &
        198, &
        199, &
        200, &
        201, &
        202, &
        203, &
        204, &
        205, &
        206, &
        207, &
        208, &
        209, &
        210, &
        211, &
        212, &
        213, &
        214, &
        215, &
        216, &
        217, &
        218, &
        219, &
        220, &
        221, &
        222, &
        223, &
        224, &
        225, &
        226, &
        227, &
        228, &
        229, &
        230, &
        231, &
        232, &
        233, &
        234, &
        235, &
        236, &
        237, &
        238, &
        239, &
        240, &
        241, &
        242, &
        243, &
        244, &
        245, &
        246, &
        247, &
        248, &
        249, &
        250, &
        251, &
        252, &
        253, &
        254, &
        255, &
        256, &
        257, &
        258, &
        259, &
        260, &
        261, &
        262, &
        263, &
        264, &
        265, &
        266, &
        267, &
        268, &
        269, &
        270, &
        271, &
        272, &
        273, &
        274, &
        275, &
        276, &
        277, &
        278, &
        279, &
        280, &
        281, &
        282, &
        283, &
        284, &
        285, &
        286, &
        287, &
        288, &
        289, &
        290, &
        291, &
        292, &
        293, &
        294, &
        295, &
        296, &
        297, &
        298, &
        299, &
        300, &
        301, &
        302, &
        303, &
        304, &
        305, &
        306, &
        307, &
        308, &
        309, &
        310, &
        311, &
        312, &
        313, &
        314, &
        315, &
        316, &
        317, &
        318, &
        319, &
        320, &
        321, &
        322, &
        323, &
        324, &
        325, &
        326, &
        327, &
        328, &
        329, &
        330, &
        331, &
        332, &
        333, &
        334, &
        335, &
        336, &
        337 &
    ]

contains

    subroutine call_loop(loop_number, in_site, out_site, A11, blocking_level, gauge_field_blocked)
        integer, intent(in) :: loop_number, in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        select case(loop_number)
        case(1)
            call loop_1(in_site, out_site, A11, blocking_level)
        case(2)
            call loop_2(in_site, out_site, A11, blocking_level)
        case(3)
            call loop_3(in_site, out_site, A11, blocking_level)
        case(4)
            call loop_4(in_site, out_site, A11, blocking_level)
        case(5)
            call loop_5(in_site, out_site, A11, blocking_level)
        case(6)
            call loop_6(in_site, out_site, A11, blocking_level)
        case(7)
            call loop_7(in_site, out_site, A11, blocking_level)
        case(8)
            call loop_8(in_site, out_site, A11, blocking_level)
        case(9)
            call loop_9(in_site, out_site, A11, blocking_level)
        case(10)
            call loop_10(in_site, out_site, A11, blocking_level)
        case(11)
            call loop_11(in_site, out_site, A11, blocking_level)
        case(12)
            call loop_12(in_site, out_site, A11, blocking_level)
        case(13)
            call loop_13(in_site, out_site, A11, blocking_level)
        case(14)
            call loop_14(in_site, out_site, A11, blocking_level)
        case(15)
            call loop_15(in_site, out_site, A11, blocking_level)
        case(16)
            call loop_16(in_site, out_site, A11, blocking_level)
        case(17)
            call loop_17(in_site, out_site, A11, blocking_level)
        case(18)
            call loop_18(in_site, out_site, A11, blocking_level)
        case(19)
            call loop_19(in_site, out_site, A11, blocking_level)
        case(20)
            call loop_20(in_site, out_site, A11, blocking_level)
        case(21)
            call loop_21(in_site, out_site, A11, blocking_level)
        case(22)
            call loop_22(in_site, out_site, A11, blocking_level)
        case(23)
            call loop_23(in_site, out_site, A11, blocking_level)
        case(24)
            call loop_24(in_site, out_site, A11, blocking_level)
        case(25)
            call loop_25(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(26)
            call loop_26(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(27)
            call loop_27(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(28)
            call loop_28(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(29)
            call loop_29(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(30)
            call loop_30(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(31)
            call loop_31(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(32)
            call loop_32(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(33)
            call loop_33(in_site, out_site, A11, blocking_level)
        case(34)
            call loop_34(in_site, out_site, A11, blocking_level)
        case(35)
            call loop_35(in_site, out_site, A11, blocking_level)
        case(36)
            call loop_36(in_site, out_site, A11, blocking_level)
        case(37)
            call loop_37(in_site, out_site, A11, blocking_level)
        case(38)
            call loop_38(in_site, out_site, A11, blocking_level)
        case(39)
            call loop_39(in_site, out_site, A11, blocking_level)
        case(40)
            call loop_40(in_site, out_site, A11, blocking_level)
        case(41)
            call loop_41(in_site, out_site, A11, blocking_level)
        case(42)
            call loop_42(in_site, out_site, A11, blocking_level)
        case(43)
            call loop_43(in_site, out_site, A11, blocking_level)
        case(44)
            call loop_44(in_site, out_site, A11, blocking_level)
        case(45)
            call loop_45(in_site, out_site, A11, blocking_level)
        case(46)
            call loop_46(in_site, out_site, A11, blocking_level)
        case(47)
            call loop_47(in_site, out_site, A11, blocking_level)
        case(48)
            call loop_48(in_site, out_site, A11, blocking_level)
        case(49)
            call loop_49(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(50)
            call loop_50(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(51)
            call loop_51(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(52)
            call loop_52(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(53)
            call loop_53(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(54)
            call loop_54(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(55)
            call loop_55(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(56)
            call loop_56(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(57)
            call loop_57(in_site, out_site, A11, blocking_level)
        case(58)
            call loop_58(in_site, out_site, A11, blocking_level)
        case(59)
            call loop_59(in_site, out_site, A11, blocking_level)
        case(60)
            call loop_60(in_site, out_site, A11, blocking_level)
        case(61)
            call loop_61(in_site, out_site, A11, blocking_level)
        case(62)
            call loop_62(in_site, out_site, A11, blocking_level)
        case(63)
            call loop_63(in_site, out_site, A11, blocking_level)
        case(64)
            call loop_64(in_site, out_site, A11, blocking_level)
        case(65)
            call loop_65(in_site, out_site, A11, blocking_level)
        case(66)
            call loop_66(in_site, out_site, A11, blocking_level)
        case(67)
            call loop_67(in_site, out_site, A11, blocking_level)
        case(68)
            call loop_68(in_site, out_site, A11, blocking_level)
        case(69)
            call loop_69(in_site, out_site, A11, blocking_level)
        case(70)
            call loop_70(in_site, out_site, A11, blocking_level)
        case(71)
            call loop_71(in_site, out_site, A11, blocking_level)
        case(72)
            call loop_72(in_site, out_site, A11, blocking_level)
        case(73)
            call loop_73(in_site, out_site, A11, blocking_level)
        case(74)
            call loop_74(in_site, out_site, A11, blocking_level)
        case(75)
            call loop_75(in_site, out_site, A11, blocking_level)
        case(76)
            call loop_76(in_site, out_site, A11, blocking_level)
        case(77)
            call loop_77(in_site, out_site, A11, blocking_level)
        case(78)
            call loop_78(in_site, out_site, A11, blocking_level)
        case(79)
            call loop_79(in_site, out_site, A11, blocking_level)
        case(80)
            call loop_80(in_site, out_site, A11, blocking_level)
        case(81)
            call loop_81(in_site, out_site, A11, blocking_level)
        case(82)
            call loop_82(in_site, out_site, A11, blocking_level)
        case(83)
            call loop_83(in_site, out_site, A11, blocking_level)
        case(84)
            call loop_84(in_site, out_site, A11, blocking_level)
        case(85)
            call loop_85(in_site, out_site, A11, blocking_level)
        case(86)
            call loop_86(in_site, out_site, A11, blocking_level)
        case(87)
            call loop_87(in_site, out_site, A11, blocking_level)
        case(88)
            call loop_88(in_site, out_site, A11, blocking_level)
        case(89)
            call loop_89(in_site, out_site, A11, blocking_level)
        case(90)
            call loop_90(in_site, out_site, A11, blocking_level)
        case(91)
            call loop_91(in_site, out_site, A11, blocking_level)
        case(92)
            call loop_92(in_site, out_site, A11, blocking_level)
        case(93)
            call loop_93(in_site, out_site, A11, blocking_level)
        case(94)
            call loop_94(in_site, out_site, A11, blocking_level)
        case(95)
            call loop_95(in_site, out_site, A11, blocking_level)
        case(96)
            call loop_96(in_site, out_site, A11, blocking_level)
        case(97)
            call loop_97(in_site, out_site, A11, blocking_level)
        case(98)
            call loop_98(in_site, out_site, A11, blocking_level)
        case(99)
            call loop_99(in_site, out_site, A11, blocking_level)
        case(100)
            call loop_100(in_site, out_site, A11, blocking_level)
        case(101)
            call loop_101(in_site, out_site, A11, blocking_level)
        case(102)
            call loop_102(in_site, out_site, A11, blocking_level)
        case(103)
            call loop_103(in_site, out_site, A11, blocking_level)
        case(104)
            call loop_104(in_site, out_site, A11, blocking_level)
        case(105)
            call loop_105(in_site, out_site, A11, blocking_level)
        case(106)
            call loop_106(in_site, out_site, A11, blocking_level)
        case(107)
            call loop_107(in_site, out_site, A11, blocking_level)
        case(108)
            call loop_108(in_site, out_site, A11, blocking_level)
        case(109)
            call loop_109(in_site, out_site, A11, blocking_level)
        case(110)
            call loop_110(in_site, out_site, A11, blocking_level)
        case(111)
            call loop_111(in_site, out_site, A11, blocking_level)
        case(112)
            call loop_112(in_site, out_site, A11, blocking_level)
        case(113)
            call loop_113(in_site, out_site, A11, blocking_level)
        case(114)
            call loop_114(in_site, out_site, A11, blocking_level)
        case(115)
            call loop_115(in_site, out_site, A11, blocking_level)
        case(116)
            call loop_116(in_site, out_site, A11, blocking_level)
        case(117)
            call loop_117(in_site, out_site, A11, blocking_level)
        case(118)
            call loop_118(in_site, out_site, A11, blocking_level)
        case(119)
            call loop_119(in_site, out_site, A11, blocking_level)
        case(120)
            call loop_120(in_site, out_site, A11, blocking_level)
        case(121)
            call loop_121(in_site, out_site, A11, blocking_level)
        case(122)
            call loop_122(in_site, out_site, A11, blocking_level)
        case(123)
            call loop_123(in_site, out_site, A11, blocking_level)
        case(124)
            call loop_124(in_site, out_site, A11, blocking_level)
        case(125)
            call loop_125(in_site, out_site, A11, blocking_level)
        case(126)
            call loop_126(in_site, out_site, A11, blocking_level)
        case(127)
            call loop_127(in_site, out_site, A11, blocking_level)
        case(128)
            call loop_128(in_site, out_site, A11, blocking_level)
        case(129)
            call loop_129(in_site, out_site, A11, blocking_level)
        case(130)
            call loop_130(in_site, out_site, A11, blocking_level)
        case(131)
            call loop_131(in_site, out_site, A11, blocking_level)
        case(132)
            call loop_132(in_site, out_site, A11, blocking_level)
        case(133)
            call loop_133(in_site, out_site, A11, blocking_level)
        case(134)
            call loop_134(in_site, out_site, A11, blocking_level)
        case(135)
            call loop_135(in_site, out_site, A11, blocking_level)
        case(136)
            call loop_136(in_site, out_site, A11, blocking_level)
        case(137)
            call loop_137(in_site, out_site, A11, blocking_level)
        case(138)
            call loop_138(in_site, out_site, A11, blocking_level)
        case(139)
            call loop_139(in_site, out_site, A11, blocking_level)
        case(140)
            call loop_140(in_site, out_site, A11, blocking_level)
        case(141)
            call loop_141(in_site, out_site, A11, blocking_level)
        case(142)
            call loop_142(in_site, out_site, A11, blocking_level)
        case(143)
            call loop_143(in_site, out_site, A11, blocking_level)
        case(144)
            call loop_144(in_site, out_site, A11, blocking_level)
        case(145)
            call loop_145(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(146)
            call loop_146(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(147)
            call loop_147(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(148)
            call loop_148(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(149)
            call loop_149(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(150)
            call loop_150(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(151)
            call loop_151(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(152)
            call loop_152(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(153)
            call loop_153(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(154)
            call loop_154(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(155)
            call loop_155(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(156)
            call loop_156(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(157)
            call loop_157(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(158)
            call loop_158(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(159)
            call loop_159(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(160)
            call loop_160(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(161)
            call loop_161(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(162)
            call loop_162(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(163)
            call loop_163(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(164)
            call loop_164(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(165)
            call loop_165(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(166)
            call loop_166(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(167)
            call loop_167(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(168)
            call loop_168(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(169)
            call loop_169(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(170)
            call loop_170(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(171)
            call loop_171(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(172)
            call loop_172(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(173)
            call loop_173(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(174)
            call loop_174(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(175)
            call loop_175(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(176)
            call loop_176(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(177)
            call loop_177(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(178)
            call loop_178(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(179)
            call loop_179(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(180)
            call loop_180(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(181)
            call loop_181(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(182)
            call loop_182(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(183)
            call loop_183(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(184)
            call loop_184(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(185)
            call loop_185(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(186)
            call loop_186(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(187)
            call loop_187(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(188)
            call loop_188(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(189)
            call loop_189(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(190)
            call loop_190(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(191)
            call loop_191(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(192)
            call loop_192(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(193)
            call loop_193(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(194)
            call loop_194(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(195)
            call loop_195(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(196)
            call loop_196(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(197)
            call loop_197(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(198)
            call loop_198(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(199)
            call loop_199(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(200)
            call loop_200(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(201)
            call loop_201(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(202)
            call loop_202(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(203)
            call loop_203(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(204)
            call loop_204(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(205)
            call loop_205(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(206)
            call loop_206(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(207)
            call loop_207(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(208)
            call loop_208(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(209)
            call loop_209(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(210)
            call loop_210(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(211)
            call loop_211(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(212)
            call loop_212(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(213)
            call loop_213(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(214)
            call loop_214(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(215)
            call loop_215(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(216)
            call loop_216(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(217)
            call loop_217(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(218)
            call loop_218(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(219)
            call loop_219(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(220)
            call loop_220(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(221)
            call loop_221(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(222)
            call loop_222(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(223)
            call loop_223(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(224)
            call loop_224(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(225)
            call loop_225(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(226)
            call loop_226(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(227)
            call loop_227(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(228)
            call loop_228(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(229)
            call loop_229(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(230)
            call loop_230(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(231)
            call loop_231(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(232)
            call loop_232(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(233)
            call loop_233(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(234)
            call loop_234(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(235)
            call loop_235(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(236)
            call loop_236(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(237)
            call loop_237(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(238)
            call loop_238(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(239)
            call loop_239(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(240)
            call loop_240(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(241)
            call loop_241(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(242)
            call loop_242(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(243)
            call loop_243(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(244)
            call loop_244(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(245)
            call loop_245(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(246)
            call loop_246(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(247)
            call loop_247(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(248)
            call loop_248(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(249)
            call loop_249(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(250)
            call loop_250(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(251)
            call loop_251(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(252)
            call loop_252(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(253)
            call loop_253(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(254)
            call loop_254(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(255)
            call loop_255(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(256)
            call loop_256(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(257)
            call loop_257(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(258)
            call loop_258(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(259)
            call loop_259(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(260)
            call loop_260(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(261)
            call loop_261(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(262)
            call loop_262(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(263)
            call loop_263(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(264)
            call loop_264(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(265)
            call loop_265(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(266)
            call loop_266(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(267)
            call loop_267(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(268)
            call loop_268(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(269)
            call loop_269(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(270)
            call loop_270(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(271)
            call loop_271(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(272)
            call loop_272(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(273)
            call loop_273(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(274)
            call loop_274(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(275)
            call loop_275(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(276)
            call loop_276(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(277)
            call loop_277(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(278)
            call loop_278(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(279)
            call loop_279(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(280)
            call loop_280(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(281)
            call loop_281(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(282)
            call loop_282(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(283)
            call loop_283(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(284)
            call loop_284(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(285)
            call loop_285(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(286)
            call loop_286(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(287)
            call loop_287(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(288)
            call loop_288(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(289)
            call loop_289(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(290)
            call loop_290(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(291)
            call loop_291(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(292)
            call loop_292(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(293)
            call loop_293(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(294)
            call loop_294(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(295)
            call loop_295(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(296)
            call loop_296(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(297)
            call loop_297(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(298)
            call loop_298(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(299)
            call loop_299(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(300)
            call loop_300(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(301)
            call loop_301(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(302)
            call loop_302(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(303)
            call loop_303(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(304)
            call loop_304(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(305)
            call loop_305(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(306)
            call loop_306(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(307)
            call loop_307(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(308)
            call loop_308(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(309)
            call loop_309(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(310)
            call loop_310(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(311)
            call loop_311(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(312)
            call loop_312(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(313)
            call loop_313(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(314)
            call loop_314(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(315)
            call loop_315(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(316)
            call loop_316(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(317)
            call loop_317(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(318)
            call loop_318(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(319)
            call loop_319(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(320)
            call loop_320(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(321)
            call loop_321(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(322)
            call loop_322(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(323)
            call loop_323(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(324)
            call loop_324(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(325)
            call loop_325(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(326)
            call loop_326(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(327)
            call loop_327(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(328)
            call loop_328(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(329)
            call loop_329(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(330)
            call loop_330(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(331)
            call loop_331(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(332)
            call loop_332(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(333)
            call loop_333(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(334)
            call loop_334(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(335)
            call loop_335(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(336)
            call loop_336(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case(337)
            call loop_337(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        case default
            error stop 'Invalid loop number in generated call_loop'
        end select
    end subroutine call_loop

end module operator_dispatch
