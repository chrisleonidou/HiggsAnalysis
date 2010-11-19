import FWCore.ParameterSet.Config as cms

tauIDFactorizationCoefficients = cms.untracked.PSet(
factorizationSourceName = cms.untracked.string('BTau_141950_148864_tctau_factorized'),

tauIDFactorizationByPt_signalAnalysisTauSelectionShrinkingConeCutBased_Coefficients = cms.untracked.vdouble( *(
  1e-99, 1e-99, 1e-99, 0.0279026, 0.017076, 0.0123371, 0.00972959, 0.00729171, 0.00512052, 0.00470804, 0.00322379, 0.00209927
) ),

tauIDFactorizationByPt_signalAnalysisTauSelectionShrinkingConeCutBased_CoefficientUncertainty = cms.untracked.vdouble( *(
  1e-99, 1e-99, 1e-99, 0.000258169, 0.000262863, 0.00029957, 0.000348599, 0.000381144, 0.000399845, 0.000367636, 0.000393848, 0.000396724
) ),

tauIDFactorizationByEta_signalAnalysisTauSelectionShrinkingConeCutBased_Coefficients = cms.untracked.vdouble( *(
  0.0195789, 0.0198828, 0.0199729, 0.0228411, 0.0269568, 0.0252487, 0.0182694, 0.0175837, 0.01681, 0.0158471, 0.014888, 0.0148763, 0.0145669, 0.0144828, 0.0157247, 0.0155326, 0.0160005, 0.0173721, 0.0165342, 0.0204219, 0.0259194, 0.0253636, 0.0203856, 0.0188299, 0.0189322, 1e-99
) ),

tauIDFactorizationByEta_signalAnalysisTauSelectionShrinkingConeCutBased_CoefficientUncertainty = cms.untracked.vdouble( *(
  0.000997834, 0.000965593, 0.000644959, 0.000739896, 0.000854159, 0.000721097, 0.000719356, 0.000963581, 0.000621741, 0.000468121, 0.000757777, 0.00051982, 0.000508698, 0.000407844, 0.000539354, 0.000774697, 0.000475985, 0.000913058, 0.000690123, 0.00080851, 0.000723904, 0.000852106, 0.000988849, 0.000542668, 0.00104535, 1e-99
) ),

tauIDFactorizationByPtVsEta_signalAnalysisTauSelectionShrinkingConeCutBased_Coefficients = cms.untracked.vdouble( *(
  1e-99, 1e-99, 1e-99, 0.0281707, 0.0179023, 0.012303, 0.0135227, 0.00632911, 0.00967742, 0.00504202, 0.00621118, 1e-99,
  1e-99, 1e-99, 1e-99, 0.0286543, 0.0191378, 0.0137615, 0.00929368, 0.00737619, 0.00497512, 0.00310559, 0.00542005, 1e-99,
  1e-99, 1e-99, 1e-99, 0.028886, 0.0186664, 0.0125816, 0.0102125, 0.0116063, 0.00585652, 0.00401875, 0.00238095, 0.00408163,
  1e-99, 1e-99, 1e-99, 0.0312934, 0.0221406, 0.0155786, 0.013306, 0.0138889, 0.00788782, 0.00408163, 0.00403226, 0.00468384,
  1e-99, 1e-99, 1e-99, 0.0398029, 0.0197986, 0.0210325, 0.0161167, 0.0104809, 0.0106383, 0.00837989, 0.00459418, 0.00520833,
  1e-99, 1e-99, 1e-99, 0.0378044, 0.022299, 0.0151351, 0.0143132, 0.0116227, 0.00639205, 0.00349406, 0.00804598, 0.00176991,
  1e-99, 1e-99, 1e-99, 0.0288523, 0.016961, 0.00926121, 0.0101302, 0.00747126, 0.000909918, 0.00580431, 0.00283688, 1e-99,
  1e-99, 1e-99, 1e-99, 0.0274428, 0.0179028, 0.0121488, 0.00625, 0.00527983, 0.00304414, 0.00149254, 0.00241546, 0.00371747,
  1e-99, 1e-99, 1e-99, 0.0273004, 0.0141948, 0.00984286, 0.00898551, 0.00528402, 0.00781805, 0.001868, 0.0010101, 0.00615385,
  1e-99, 1e-99, 1e-99, 0.0249481, 0.015056, 0.00961733, 0.00862218, 0.00536193, 0.00485633, 0.00376648, 0.00183374, 0.000931966,
  1e-99, 1e-99, 1e-99, 0.0219186, 0.0152259, 0.0107072, 0.00813787, 0.003125, 0.00367197, 0.0052687, 0.00504202, 1e-99,
  1e-99, 1e-99, 1e-99, 0.0226514, 0.0140781, 0.0119672, 0.00800712, 0.003861, 0.00334635, 0.00254194, 0.00307692, 0.00126263,
  1e-99, 1e-99, 1e-99, 0.023604, 0.0126659, 0.00931048, 0.00652032, 0.00449205, 0.0021846, 0.00449775, 0.000801925, 0.00124844,
  1e-99, 1e-99, 1e-99, 0.0228909, 0.0131875, 0.00997801, 0.00743069, 0.00555556, 0.00105411, 0.00343643, 0.00161812, 1e-99,
  1e-99, 1e-99, 1e-99, 0.0244597, 0.015387, 0.0105451, 0.00664985, 0.00333704, 0.00545852, 0.00351229, 0.00332502, 0.00343249,
  1e-99, 1e-99, 1e-99, 0.024777, 0.0153846, 0.0111429, 0.00471254, 0.00521221, 1e-99, 0.00106838, 0.00351494, 1e-99,
  1e-99, 1e-99, 1e-99, 0.0241158, 0.0157286, 0.0113188, 0.00915493, 0.00573301, 0.00521966, 0.00333087, 0.00324675, 0.00087184,
  1e-99, 1e-99, 1e-99, 0.0287462, 0.0146586, 0.0100466, 0.00622665, 0.00657277, 0.00421348, 0.00662252, 0.00219298, 1e-99,
  1e-99, 1e-99, 1e-99, 0.0274929, 0.0133803, 0.00962567, 0.00664943, 0.00782123, 0.00183655, 0.00640512, 0.00136986, 0.00401606,
  1e-99, 1e-99, 1e-99, 0.0292388, 0.0203363, 0.0155794, 0.0108878, 0.00526662, 0.00727651, 0.00477555, 0.00505902, 0.0027248,
  1e-99, 1e-99, 1e-99, 0.0376476, 0.0224129, 0.0165789, 0.0154278, 0.0126464, 0.00668151, 0.00775194, 0.00464037, 0.00591716,
  1e-99, 1e-99, 1e-99, 0.033217, 0.0245196, 0.0184502, 0.0156062, 0.0137363, 0.00797267, 0.0149733, 0.00743494, 1e-99,
  1e-99, 1e-99, 1e-99, 0.0280272, 0.0185644, 0.0142805, 0.0119681, 0.0124717, 0.00883392, 0.00661157, 0.00265252, 0.00497512,
  1e-99, 1e-99, 1e-99, 0.026464, 0.0179234, 0.0122268, 0.00985532, 0.00927198, 0.00802139, 0.00821355, 0.00360036, 0.00316456,
  1e-99, 1e-99, 1e-99, 0.0248838, 0.0186486, 0.0157411, 0.0126771, 0.00739827, 0.0117647, 0.00568182, 0.00377358, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99
) ),

tauIDFactorizationByPtVsEta_signalAnalysisTauSelectionShrinkingConeCutBased_CoefficientUncertainty = cms.untracked.vdouble( *(
  1e-99, 1e-99, 1e-99, 0.0018535, 0.00194177, 0.00217488, 0.00302376, 0.00258385, 0.00395079, 0.00291101, 0.00439197, 0.00662252,
  1e-99, 1e-99, 1e-99, 0.0017874, 0.00192342, 0.0022036, 0.00239962, 0.00278794, 0.00287239, 0.00219598, 0.00383256, 0.00588235,
  1e-99, 1e-99, 1e-99, 0.00119429, 0.00127304, 0.00141554, 0.00167893, 0.00232126, 0.00207059, 0.00164065, 0.00168359, 0.00288615,
  1e-99, 1e-99, 1e-99, 0.00131886, 0.00148264, 0.00169977, 0.00213067, 0.00277778, 0.00262927, 0.00182536, 0.00232803, 0.00331198,
  1e-99, 1e-99, 1e-99, 0.00157581, 0.00149663, 0.00211385, 0.00248685, 0.00254199, 0.00320757, 0.0027933, 0.00265245, 0.00368285,
  1e-99, 1e-99, 1e-99, 0.00136238, 0.00138828, 0.00152888, 0.00198489, 0.0022794, 0.00213068, 0.00156259, 0.00304109, 0.00176991,
  1e-99, 1e-99, 1e-99, 0.00143368, 0.00140854, 0.00139618, 0.00191444, 0.00207216, 0.000909918, 0.00219382, 0.00200598, 0.00220264,
  1e-99, 1e-99, 1e-99, 0.00195028, 0.00195336, 0.00214763, 0.00208333, 0.00236121, 0.00215253, 0.00149254, 0.00241546, 0.00371747,
  1e-99, 1e-99, 1e-99, 0.00126331, 0.00117881, 0.00130372, 0.00161384, 0.00152536, 0.00235723, 0.00107849, 0.0010101, 0.00307692,
  1e-99, 1e-99, 1e-99, 0.000944977, 0.000937347, 0.000986717, 0.00121936, 0.00119896, 0.0014019, 0.00119107, 0.00105871, 0.000931966,
  1e-99, 1e-99, 1e-99, 0.00146778, 0.00157886, 0.00173694, 0.00197372, 0.0015625, 0.00212001, 0.00235624, 0.00291101, 0.00273224,
  1e-99, 1e-99, 1e-99, 0.00103174, 0.00103504, 0.00126852, 0.00133452, 0.00116414, 0.00136614, 0.00113679, 0.00153846, 0.00126263,
  1e-99, 1e-99, 1e-99, 0.00103016, 0.000977193, 0.00112085, 0.00119044, 0.00124587, 0.0010923, 0.00149925, 0.000801925, 0.00124844,
  1e-99, 1e-99, 1e-99, 0.000820151, 0.000802563, 0.000918551, 0.00103045, 0.00111111, 0.000608591, 0.00103612, 0.000934224, 0.000746826,
  1e-99, 1e-99, 1e-99, 0.00107783, 0.00109907, 0.00120172, 0.00123485, 0.00111235, 0.00172613, 0.00132752, 0.00166251, 0.00198175,
  1e-99, 1e-99, 1e-99, 0.00156704, 0.00159531, 0.00178429, 0.00149023, 0.00197003, 0.00114943, 0.00106838, 0.00248544, 0.00246305,
  1e-99, 1e-99, 1e-99, 0.000938705, 0.000973574, 0.00108414, 0.00126956, 0.00125104, 0.00150679, 0.00111029, 0.00145199, 0.00087184,
  1e-99, 1e-99, 1e-99, 0.00187519, 0.00171566, 0.00189864, 0.00196904, 0.00248427, 0.00243266, 0.00296168, 0.00219298, 0.00331126,
  1e-99, 1e-99, 1e-99, 0.00142929, 0.00125318, 0.00143491, 0.00156729, 0.00209031, 0.00129864, 0.00226455, 0.00136986, 0.00283979,
  1e-99, 1e-99, 1e-99, 0.00152418, 0.00162821, 0.00194742, 0.00213527, 0.00186203, 0.00275026, 0.00213569, 0.00292083, 0.0027248,
  1e-99, 1e-99, 1e-99, 0.00133355, 0.00136654, 0.00160274, 0.00208028, 0.0024338, 0.00222717, 0.0023373, 0.00232019, 0.00341627,
  1e-99, 1e-99, 1e-99, 0.00146514, 0.00170013, 0.00206279, 0.002499, 0.00307152, 0.00301338, 0.00400177, 0.00371747, 0.003003,
  1e-99, 1e-99, 1e-99, 0.0017517, 0.00195685, 0.00228671, 0.0028209, 0.00376035, 0.00395065, 0.00330579, 0.00265252, 0.00497512,
  1e-99, 1e-99, 1e-99, 0.000990392, 0.00107113, 0.00122883, 0.00143755, 0.00178439, 0.00207111, 0.00205339, 0.00180018, 0.00223768,
  1e-99, 1e-99, 1e-99, 0.00184451, 0.0021252, 0.00262352, 0.00307465, 0.00302033, 0.00480292, 0.0032804, 0.00377358, 0.00714286,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99
) ),

tauIDFactorizationByPt_signalAnalysisTauSelectionShrinkingConeTaNCBased_Coefficients = cms.untracked.vdouble( *(
  1e-99, 1e-99, 1e-99, 0.0184074, 0.0109092, 0.00797253, 0.00608256, 0.0050006, 0.00346572, 0.0034162, 0.00245393, 0.00134953
) ),

tauIDFactorizationByPt_signalAnalysisTauSelectionShrinkingConeTaNCBased_CoefficientUncertainty = cms.untracked.vdouble( *(
  1e-99, 1e-99, 1e-99, 0.000209691, 0.000210103, 0.000240819, 0.000275627, 0.000315635, 0.000328951, 0.000313163, 0.000343619, 0.000318087
) ),

tauIDFactorizationByEta_signalAnalysisTauSelectionShrinkingConeTaNCBased_Coefficients = cms.untracked.vdouble( *(
  0.00635679, 0.00886284, 0.0104134, 0.0139731, 0.0172404, 0.0170933, 0.0128027, 0.0127257, 0.0118659, 0.0113944, 0.011031, 0.0111164, 0.0104455, 0.010417, 0.0114143, 0.011437, 0.0116109, 0.0117094, 0.0115221, 0.013796, 0.0176301, 0.0156017, 0.01132, 0.0102282, 0.00750361, 1e-99
) ),

tauIDFactorizationByEta_signalAnalysisTauSelectionShrinkingConeTaNCBased_CoefficientUncertainty = cms.untracked.vdouble( *(
  0.000568569, 0.000644677, 0.000465702, 0.000578707, 0.000683091, 0.000593318, 0.000602189, 0.000819737, 0.000522367, 0.000396944, 0.000652275, 0.000449352, 0.000430767, 0.000345892, 0.000459522, 0.00066476, 0.000405472, 0.000749616, 0.000576103, 0.000664529, 0.000597029, 0.000668305, 0.000736871, 0.000399953, 0.00065811, 1e-99
) ),

tauIDFactorizationByPtVsEta_signalAnalysisTauSelectionShrinkingConeTaNCBased_Coefficients = cms.untracked.vdouble( *(
  1e-99, 1e-99, 1e-99, 0.00756098, 0.00652906, 0.00499808, 0.00540906, 0.00316456, 0.00322581, 0.00840336, 0.00310559, 1e-99,
  1e-99, 1e-99, 1e-99, 0.011707, 0.00908564, 0.00635145, 0.00495663, 0.00737619, 0.00165837, 0.00310559, 0.00271003, 1e-99,
  1e-99, 1e-99, 1e-99, 0.015159, 0.00972391, 0.00668896, 0.00358819, 0.00649954, 0.00366032, 0.00334896, 0.00119048, 0.00204082,
  1e-99, 1e-99, 1e-99, 0.0187316, 0.0131056, 0.0100148, 0.00784715, 0.00944444, 0.00876424, 0.00408163, 0.00537634, 0.00234192,
  1e-99, 1e-99, 1e-99, 0.0254539, 0.012558, 0.0138092, 0.00997698, 0.00739827, 0.00676983, 0.00465549, 0.00153139, 0.00520833,
  1e-99, 1e-99, 1e-99, 0.0263158, 0.014261, 0.00926641, 0.00990917, 0.00849352, 0.00426136, 0.00139762, 0.00689655, 1e-99,
  1e-99, 1e-99, 1e-99, 0.0207309, 0.0111124, 0.00610398, 0.0068741, 0.00402299, 0.00181984, 0.0066335, 0.00141844, 1e-99,
  1e-99, 1e-99, 1e-99, 0.0196812, 0.0125746, 0.00911162, 0.00694444, 0.0031679, 0.00152207, 1e-99, 0.00483092, 1e-99,
  1e-99, 1e-99, 1e-99, 0.0189992, 0.0109643, 0.00587118, 0.00608696, 0.00484368, 0.00497512, 0.001868, 0.0030303, 1e-99,
  1e-99, 1e-99, 1e-99, 0.0185768, 0.0102708, 0.00647904, 0.00551819, 0.00268097, 0.00404694, 0.00338983, 0.00183374, 0.000931966,
  1e-99, 1e-99, 1e-99, 0.0159229, 0.0111329, 0.00873485, 0.00622307, 0.00546875, 1e-99, 0.00421496, 0.00168067, 1e-99,
  1e-99, 1e-99, 1e-99, 0.0172001, 0.0105776, 0.00820223, 0.00511566, 0.002808, 0.0044618, 0.00203355, 0.00153846, 0.00126263,
  1e-99, 1e-99, 1e-99, 0.0172197, 0.00859469, 0.00769127, 0.0034775, 0.00310988, 0.0021846, 0.00149925, 0.000801925, 0.00124844,
  1e-99, 1e-99, 1e-99, 0.0168376, 0.00844974, 0.00735667, 0.00585882, 0.00488889, 0.00035137, 0.00187441, 0.0021575, 1e-99,
  1e-99, 1e-99, 1e-99, 0.0178105, 0.0114618, 0.00739523, 0.0045861, 0.00259548, 0.00327511, 0.00250878, 0.00166251, 0.00228833,
  1e-99, 1e-99, 1e-99, 0.018335, 0.00992556, 0.00914286, 0.00471254, 0.00446761, 1e-99, 0.00213675, 0.00175747, 1e-99,
  1e-99, 1e-99, 1e-99, 0.0175387, 0.0121731, 0.00809969, 0.00580986, 0.003822, 0.00173989, 0.00185048, 0.00194805, 0.00087184,
  1e-99, 1e-99, 1e-99, 0.0192049, 0.0108434, 0.00574094, 0.00373599, 0.0028169, 0.00421348, 0.00529801, 1e-99, 0.00331126,
  1e-99, 1e-99, 1e-99, 0.0199138, 0.00833333, 0.0059893, 0.00443295, 0.00670391, 0.00183655, 0.00320256, 0.00273973, 0.00200803,
  1e-99, 1e-99, 1e-99, 0.0210551, 0.0132968, 0.00851996, 0.00544389, 0.00263331, 0.00519751, 0.00286533, 0.00337268, 0.00544959,
  1e-99, 1e-99, 1e-99, 0.0266887, 0.0144143, 0.0114658, 0.00813464, 0.00655738, 0.00445434, 0.00493305, 0.00348028, 0.00197239,
  1e-99, 1e-99, 1e-99, 0.0210676, 0.0140281, 0.0103782, 0.0104042, 0.00824176, 0.00569476, 0.00962567, 0.00371747, 0.003003,
  1e-99, 1e-99, 1e-99, 0.0153273, 0.00990099, 0.00768949, 0.00731383, 0.00793651, 0.00353357, 0.00826446, 0.00265252, 0.00497512,
  1e-99, 1e-99, 1e-99, 0.0138992, 0.00972987, 0.00716315, 0.00566156, 0.00583791, 0.00534759, 0.00564682, 0.00270027, 0.00158228,
  1e-99, 1e-99, 1e-99, 0.0073831, 0.00847663, 0.00699606, 0.00820283, 0.00739827, 0.00784314, 0.00568182, 0.00377358, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99
) ),

tauIDFactorizationByPtVsEta_signalAnalysisTauSelectionShrinkingConeTaNCBased_CoefficientUncertainty = cms.untracked.vdouble( *(
  1e-99, 1e-99, 1e-99, 0.000960245, 0.00117265, 0.00138622, 0.00191239, 0.00182706, 0.00228099, 0.0037581, 0.00310559, 0.00662252,
  1e-99, 1e-99, 1e-99, 0.00114249, 0.00132528, 0.00149705, 0.00175243, 0.00278794, 0.00165837, 0.00219598, 0.00271003, 0.00588235,
  1e-99, 1e-99, 1e-99, 0.00086517, 0.000918823, 0.00103213, 0.000995184, 0.00173707, 0.00163695, 0.0014977, 0.00119048, 0.00204082,
  1e-99, 1e-99, 1e-99, 0.00102037, 0.0011407, 0.00136285, 0.00163624, 0.00229061, 0.0027715, 0.00182536, 0.00268817, 0.00234192,
  1e-99, 1e-99, 1e-99, 0.00126015, 0.00119195, 0.00171282, 0.00195665, 0.0021357, 0.00255875, 0.002082, 0.00153139, 0.00368285,
  1e-99, 1e-99, 1e-99, 0.00113667, 0.00111022, 0.00119629, 0.00165153, 0.00194855, 0.00173969, 0.000988269, 0.00281551, 0.00176991,
  1e-99, 1e-99, 1e-99, 0.00121527, 0.00114011, 0.00113348, 0.00157703, 0.00152055, 0.00128682, 0.0023453, 0.00141844, 0.00220264,
  1e-99, 1e-99, 1e-99, 0.00165161, 0.00163707, 0.0018599, 0.00219603, 0.00182899, 0.00152207, 0.00149254, 0.00341597, 0.00371747,
  1e-99, 1e-99, 1e-99, 0.00105388, 0.00103603, 0.0010069, 0.00132828, 0.00146042, 0.00188042, 0.00107849, 0.00174955, 0.00153846,
  1e-99, 1e-99, 1e-99, 0.000815433, 0.000774189, 0.000809881, 0.000975488, 0.000847796, 0.00127976, 0.00112994, 0.00105871, 0.000931966,
  1e-99, 1e-99, 1e-99, 0.00125102, 0.00135007, 0.00156883, 0.00172597, 0.00206699, 0.00122399, 0.00210748, 0.00168067, 0.00273224,
  1e-99, 1e-99, 1e-99, 0.000899061, 0.000897179, 0.00105019, 0.00106669, 0.000992779, 0.00157748, 0.00101678, 0.00108786, 0.00126263,
  1e-99, 1e-99, 1e-99, 0.000879884, 0.000804967, 0.00101873, 0.000869376, 0.00103663, 0.0010923, 0.000865593, 0.000801925, 0.00124844,
  1e-99, 1e-99, 1e-99, 0.0007034, 0.000642422, 0.000788718, 0.000914993, 0.00104231, 0.00035137, 0.000765226, 0.00107875, 0.000746826,
  1e-99, 1e-99, 1e-99, 0.00091973, 0.000948583, 0.00100636, 0.00102548, 0.000980998, 0.00133706, 0.00112196, 0.00117557, 0.00161809,
  1e-99, 1e-99, 1e-99, 0.00134801, 0.00128138, 0.00161624, 0.00149023, 0.00182389, 0.00114943, 0.00151091, 0.00175747, 0.00246305,
  1e-99, 1e-99, 1e-99, 0.00080053, 0.000856495, 0.000917109, 0.00101137, 0.00102147, 0.000869943, 0.00082756, 0.00112471, 0.00087184,
  1e-99, 1e-99, 1e-99, 0.00153272, 0.0014756, 0.00143524, 0.00152521, 0.00162634, 0.00243266, 0.00264901, 0.00219298, 0.00331126,
  1e-99, 1e-99, 1e-99, 0.00121643, 0.000988985, 0.00113187, 0.00127968, 0.00193525, 0.00129864, 0.00160128, 0.00193728, 0.00200803,
  1e-99, 1e-99, 1e-99, 0.00129341, 0.00131658, 0.00144014, 0.00150986, 0.00131666, 0.00232439, 0.0016543, 0.00238485, 0.00385344,
  1e-99, 1e-99, 1e-99, 0.0011228, 0.0010959, 0.00133287, 0.00151057, 0.00175253, 0.00181848, 0.00186452, 0.00200934, 0.00197239,
  1e-99, 1e-99, 1e-99, 0.00116683, 0.00128595, 0.0015471, 0.00204042, 0.00237919, 0.00254677, 0.00320856, 0.00262865, 0.003003,
  1e-99, 1e-99, 1e-99, 0.0012954, 0.00142908, 0.00167798, 0.0022052, 0.00299972, 0.00249861, 0.00369598, 0.00265252, 0.00497512,
  1e-99, 1e-99, 1e-99, 0.000717751, 0.000789197, 0.000940567, 0.00108957, 0.0014159, 0.00169106, 0.00170258, 0.001559, 0.00158228,
  1e-99, 1e-99, 1e-99, 0.00100471, 0.00143281, 0.00174902, 0.00247325, 0.00302033, 0.00392157, 0.0032804, 0.00377358, 0.00714286,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99
) ),

tauIDFactorizationByPt_signalAnalysisTauSelectionHPSTauBased_Coefficients = cms.untracked.vdouble( *(
  1e-99, 1e-99, 1e-99, 0.00903636, 0.00602056, 0.00413281, 0.00365002, 0.00282591, 0.00162551, 0.0022089, 0.00191514, 0.000553557
) ),

tauIDFactorizationByPt_signalAnalysisTauSelectionHPSTauBased_CoefficientUncertainty = cms.untracked.vdouble( *(
  1e-99, 1e-99, 1e-99, 0.000177286, 0.000208851, 0.000241031, 0.000309591, 0.00036182, 0.000354716, 0.000417443, 0.000531165, 0.000391424
) ),

tauIDFactorizationByEta_signalAnalysisTauSelectionHPSTauBased_Coefficients = cms.untracked.vdouble( *(
  0.00687081, 0.00948648, 0.00864492, 0.0109097, 0.0117833, 0.00978932, 0.00515009, 0.00493376, 0.00484962, 0.00472321, 0.00491383, 0.00454796, 0.004312, 0.00415524, 0.00488901, 0.0053719, 0.004243, 0.00529012, 0.0049562, 0.00649279, 0.0107698, 0.0112695, 0.00972955, 0.00867466, 0.00714769, 1e-99
) ),

tauIDFactorizationByEta_signalAnalysisTauSelectionHPSTauBased_CoefficientUncertainty = cms.untracked.vdouble( *(
  0.000763424, 0.000862407, 0.000554572, 0.000672721, 0.000739348, 0.000589248, 0.000502597, 0.0006714, 0.000444564, 0.000339984, 0.000583165, 0.000384373, 0.000371118, 0.000291641, 0.000401874, 0.000608248, 0.000325423, 0.000666492, 0.000498116, 0.000600258, 0.000609718, 0.000754659, 0.000895678, 0.00047971, 0.000830902, 1e-99
) ),

tauIDFactorizationByPtVsEta_signalAnalysisTauSelectionHPSTauBased_Coefficients = cms.untracked.vdouble( *(
  1e-99, 1e-99, 1e-99, 0.0101288, 0.00547445, 0.00136519, 0.00131752, 0.0075188, 1e-99, 0.00471698, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 0.0116618, 0.00979752, 0.00805452, 0.00245098, 0.00447427, 0.00387597, 1e-99, 0.00763359, 1e-99,
  1e-99, 1e-99, 1e-99, 0.0110116, 0.00930576, 0.00388292, 0.00458716, 0.00487329, 0.00323625, 0.00172117, 1e-99, 0.00840336,
  1e-99, 1e-99, 1e-99, 0.0144682, 0.00883392, 0.00802512, 0.00579151, 0.00969697, 0.00204499, 0.00201207, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 0.01471, 0.0103469, 0.00811751, 0.00990854, 0.00508906, 0.00894855, 0.00647948, 0.00404858, 1e-99,
  1e-99, 1e-99, 1e-99, 0.0140298, 0.00686499, 0.00482509, 0.00723831, 0.00405268, 0.00170648, 0.0017452, 0.00303951, 1e-99,
  1e-99, 1e-99, 1e-99, 0.00721298, 0.00494943, 0.00244002, 0.00235664, 0.00130548, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 0.00694043, 0.00475059, 0.00362845, 0.00142248, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 0.00715929, 0.00418556, 0.00244073, 0.00246305, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 0.00646376, 0.00411609, 0.00304198, 0.00226672, 0.00197368, 0.00108342, 0.00105485, 0.00204918, 1e-99,
  1e-99, 1e-99, 1e-99, 0.00673642, 0.00417163, 0.00412493, 0.00204499, 1e-99, 1e-99, 1e-99, 0.00546448, 1e-99,
  1e-99, 1e-99, 1e-99, 0.00622673, 0.00320513, 0.00291777, 0.00417101, 0.000848896, 0.0027933, 0.00149477, 0.00273224, 0.00465116,
  1e-99, 1e-99, 1e-99, 0.00557231, 0.00413337, 0.00304709, 0.00284765, 0.000898473, 1e-99, 0.00299401, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 0.00576491, 0.0035233, 0.00283949, 0.00191632, 0.00168161, 1e-99, 0.000912409, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 0.00655896, 0.00479549, 0.00300958, 0.00353893, 1e-99, 1e-99, 1e-99, 0.00257069, 1e-99,
  1e-99, 1e-99, 1e-99, 0.00706561, 0.00650888, 0.00283768, 1e-99, 0.00183824, 1e-99, 0.00301205, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 0.00562418, 0.004079, 0.00304321, 0.00264051, 0.000682128, 0.00103413, 1e-99, 0.00195312, 1e-99,
  1e-99, 1e-99, 1e-99, 0.00677083, 0.00436681, 0.00363901, 0.00255102, 0.00639659, 0.00367647, 0.00367647, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 0.00755192, 0.00368364, 0.001638, 0.000784314, 0.00395778, 0.00452489, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 0.00829734, 0.00585652, 0.00559701, 0.00432152, 1e-99, 0.00245098, 0.00251889, 0.0045045, 1e-99,
  1e-99, 1e-99, 1e-99, 0.0144437, 0.00975244, 0.00468796, 0.00552181, 0.00603015, 0.00337838, 0.0034662, 0.00671141, 1e-99,
  1e-99, 1e-99, 1e-99, 0.0151116, 0.00958815, 0.00687581, 0.00473934, 0.00291121, 0.00262467, 0.0104439, 0.00505051, 1e-99,
  1e-99, 1e-99, 1e-99, 0.0129086, 0.00836528, 0.00624566, 0.00530504, 0.00239808, 1e-99, 0.0127119, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 0.0111589, 0.00817561, 0.0055054, 0.00372517, 0.00651702, 0.00123305, 0.00554785, 0.00274725, 1e-99,
  1e-99, 1e-99, 1e-99, 0.00777605, 0.00763973, 0.00741351, 0.00934579, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99
) ),

tauIDFactorizationByPtVsEta_signalAnalysisTauSelectionHPSTauBased_CoefficientUncertainty = cms.untracked.vdouble( *(
  1e-99, 1e-99, 1e-99, 0.00131865, 0.0014135, 0.000965333, 0.00131752, 0.00434098, 0.00438596, 0.00471698, 0.00877193, 0.0212766,
  1e-99, 1e-99, 1e-99, 0.00137436, 0.00178877, 0.00223392, 0.0017331, 0.00316379, 0.00387597, 0.00478469, 0.00763359, 0.0227273,
  1e-99, 1e-99, 1e-99, 0.000899093, 0.00117242, 0.00107693, 0.0016218, 0.0021794, 0.00228837, 0.00172117, 0.00355872, 0.00840336,
  1e-99, 1e-99, 1e-99, 0.00110641, 0.00124931, 0.00167335, 0.0019305, 0.0034284, 0.00204499, 0.00201207, 0.00373134, 0.00775194,
  1e-99, 1e-99, 1e-99, 0.00117399, 0.00144886, 0.00177139, 0.00274813, 0.00254453, 0.00447427, 0.00374093, 0.00404858, 0.00892857,
  1e-99, 1e-99, 1e-99, 0.00100469, 0.00102337, 0.00120627, 0.00200755, 0.00202634, 0.00170648, 0.0017452, 0.00303951, 0.00653595,
  1e-99, 1e-99, 1e-99, 0.000850058, 0.00103203, 0.000996132, 0.00136061, 0.00130548, 0.00233645, 0.00224215, 0.00413223, 0.00689655,
  1e-99, 1e-99, 1e-99, 0.00115674, 0.00137138, 0.00162269, 0.00142248, 0.00235294, 0.00392157, 0.0041841, 0.00724638, 0.0106383,
  1e-99, 1e-99, 1e-99, 0.000781143, 0.000854374, 0.000922507, 0.00123153, 0.00106724, 0.00182482, 0.00163132, 0.00321543, 0.00588235,
  1e-99, 1e-99, 1e-99, 0.000573566, 0.000659103, 0.000785436, 0.000925383, 0.00113951, 0.00108342, 0.00105485, 0.00204918, 0.0035461,
  1e-99, 1e-99, 1e-99, 0.000982608, 0.00111492, 0.00155908, 0.00144603, 0.00192678, 0.0030303, 0.00321543, 0.00546448, 0.0102041,
  1e-99, 1e-99, 1e-99, 0.000649182, 0.000668315, 0.000879741, 0.00147468, 0.000848896, 0.00197516, 0.00149477, 0.00273224, 0.00465116,
  1e-99, 1e-99, 1e-99, 0.000604402, 0.000754647, 0.000918733, 0.00116255, 0.000898473, 0.0013947, 0.00211709, 0.00277778, 0.00452489,
  1e-99, 1e-99, 1e-99, 0.000494337, 0.000557082, 0.000688676, 0.000782335, 0.00097088, 0.000993049, 0.000912409, 0.00176678, 0.00294985,
  1e-99, 1e-99, 1e-99, 0.000672935, 0.000822419, 0.000907421, 0.00133759, 0.000903342, 0.00144509, 0.00151515, 0.00257069, 0.00460829,
  1e-99, 1e-99, 1e-99, 0.00100937, 0.0013877, 0.00126905, 0.0010142, 0.00183824, 0.00316456, 0.00301205, 0.00617284, 0.00970874,
  1e-99, 1e-99, 1e-99, 0.00054371, 0.000661702, 0.000785754, 0.00099802, 0.000682128, 0.00103413, 0.0010929, 0.00195312, 0.00350877,
  1e-99, 1e-99, 1e-99, 0.0010842, 0.00126059, 0.00162741, 0.00180384, 0.00369307, 0.00367647, 0.00367647, 0.00729927, 0.0107527,
  1e-99, 1e-99, 1e-99, 0.000890002, 0.000893414, 0.000819001, 0.000784314, 0.00228503, 0.00319958, 0.00203252, 0.003663, 0.00694444,
  1e-99, 1e-99, 1e-99, 0.00097113, 0.00119546, 0.00161572, 0.00193264, 0.00146413, 0.00245098, 0.00251889, 0.0045045, 0.00884956,
  1e-99, 1e-99, 1e-99, 0.00099909, 0.00120964, 0.00117199, 0.00174615, 0.0024618, 0.00238887, 0.00245098, 0.00474568, 0.00671141,
  1e-99, 1e-99, 1e-99, 0.00123799, 0.00144547, 0.00171895, 0.00193483, 0.00205854, 0.00262467, 0.00522193, 0.00505051, 0.0103093,
  1e-99, 1e-99, 1e-99, 0.00147108, 0.00170756, 0.00208189, 0.00265252, 0.00239808, 0.00369004, 0.0073392, 0.00787402, 0.0208333,
  1e-99, 1e-99, 1e-99, 0.000779371, 0.000956882, 0.00110108, 0.00124172, 0.00217234, 0.00123305, 0.00277393, 0.00274725, 0.00617284,
  1e-99, 1e-99, 1e-99, 0.0012295, 0.00175267, 0.00247117, 0.0038154, 0.0027248, 0.00460829, 0.0060241, 0.0121951, 0.0294118,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99
) ),

tauIDFactorizationByPt_signalAnalysisTauSelectionCaloTauCutBased_Coefficients = cms.untracked.vdouble( *(
  1e-99, 1e-99, 1e-99, 0.0357284, 0.0213645, 0.0156083, 0.0118201, 0.0097621, 0.00885631, 0.00766317, 0.00667622, 0.0109338
) ),

tauIDFactorizationByPt_signalAnalysisTauSelectionCaloTauCutBased_CoefficientUncertainty = cms.untracked.vdouble( *(
  1e-99, 1e-99, 1e-99, 0.000301016, 0.000284204, 0.000300884, 0.000324479, 0.000361311, 0.000421729, 0.000364913, 0.000424798, 0.000637674
) ),

tauIDFactorizationByEta_signalAnalysisTauSelectionCaloTauCutBased_Coefficients = cms.untracked.vdouble( *(
  0.0267652, 0.0263873, 0.0276941, 0.0278868, 0.0315302, 0.018303, 0.020912, 0.0229766, 0.0190166, 0.0191979, 0.0141463, 0.019262, 0.0174685, 0.0181804, 0.0164941, 0.018924, 0.0185984, 0.0225784, 0.0189698, 0.025872, 0.0205569, 0.0308865, 0.0272261, 0.0263046, 0.0280673, 1e-99
) ),

tauIDFactorizationByEta_signalAnalysisTauSelectionCaloTauCutBased_CoefficientUncertainty = cms.untracked.vdouble( *(
  0.00109819, 0.00102713, 0.00070871, 0.000747446, 0.000839687, 0.000567278, 0.000723257, 0.0010401, 0.000624586, 0.000485751, 0.000693579, 0.000553514, 0.000519428, 0.000426037, 0.000515694, 0.000804004, 0.000482789, 0.000986343, 0.00070019, 0.000861442, 0.00058543, 0.000856307, 0.00105262, 0.000599383, 0.00120337, 1e-99
) ),

tauIDFactorizationByPtVsEta_signalAnalysisTauSelectionCaloTauCutBased_Coefficients = cms.untracked.vdouble( *(
  1e-99, 1e-99, 1e-99, 0.0399704, 0.0250988, 0.0191403, 0.0141602, 0.0101642, 0.0144578, 0.0105758, 0.0213592, 0.0252366,
  1e-99, 1e-99, 1e-99, 0.0387363, 0.026437, 0.0179021, 0.0174902, 0.010043, 0.00864865, 0.0110442, 0.00803859, 0.0408163,
  1e-99, 1e-99, 1e-99, 0.0421393, 0.0258429, 0.0187365, 0.0143354, 0.0143347, 0.0120664, 0.0145388, 0.0096225, 0.0236559,
  1e-99, 1e-99, 1e-99, 0.0445114, 0.0262242, 0.0187415, 0.0166163, 0.0145119, 0.0108055, 0.00964489, 0.0121131, 0.024777,
  1e-99, 1e-99, 1e-99, 0.0515897, 0.0283227, 0.0240148, 0.014652, 0.0143215, 0.0183844, 0.0139147, 0.0122655, 0.0169492,
  1e-99, 1e-99, 1e-99, 0.028486, 0.0194851, 0.0131499, 0.0109426, 0.0083815, 0.00562284, 0.00599475, 0.0079096, 0.0085592,
  1e-99, 1e-99, 1e-99, 0.0343036, 0.0182519, 0.0150969, 0.0125668, 0.0126034, 0.00925926, 0.0105634, 0.00626468, 0.0144928,
  1e-99, 1e-99, 1e-99, 0.0387097, 0.0209613, 0.0183196, 0.0101368, 0.0119581, 0.0154959, 0.0104803, 0.00144509, 0.0158172,
  1e-99, 1e-99, 1e-99, 0.0328641, 0.018106, 0.0135772, 0.00909668, 0.0100946, 0.0103627, 0.00359281, 0.00739477, 0.00493016,
  1e-99, 1e-99, 1e-99, 0.0337927, 0.0180373, 0.0135602, 0.0121173, 0.00766498, 0.00725446, 0.00460511, 0.00497159, 0.00605144,
  1e-99, 1e-99, 1e-99, 0.0235432, 0.0158403, 0.00901734, 0.00910364, 0.00697906, 0.00308642, 0.00395778, 0.00512821, 0.00136426,
  1e-99, 1e-99, 1e-99, 0.0333451, 0.0199612, 0.0141944, 0.00899743, 0.00773694, 0.00559701, 0.00520993, 0.00668258, 0.00468855,
  1e-99, 1e-99, 1e-99, 0.0314523, 0.0171785, 0.0121574, 0.00775685, 0.00532161, 0.00676879, 0.00440806, 0.00234522, 0.00401338,
  1e-99, 1e-99, 1e-99, 0.0325178, 0.0179859, 0.0115244, 0.0084444, 0.00885767, 0.00478142, 0.00565082, 0.00271657, 0.00587002,
  1e-99, 1e-99, 1e-99, 0.0281747, 0.0177902, 0.0111787, 0.00782811, 0.00625626, 0.00578871, 0.00615006, 0.00245218, 0.00590551,
  1e-99, 1e-99, 1e-99, 0.0344012, 0.0187003, 0.0142857, 0.00701107, 0.00581703, 0.00442804, 0.00307125, 0.00420168, 0.0106667,
  1e-99, 1e-99, 1e-99, 0.0311634, 0.0192174, 0.0133685, 0.0113548, 0.00892164, 0.00660975, 0.00553018, 0.00432432, 0.0073278,
  1e-99, 1e-99, 1e-99, 0.0393364, 0.0201839, 0.0167845, 0.0136488, 0.0125, 0.0112254, 0.00502092, 0.0025974, 0.00515464,
  1e-99, 1e-99, 1e-99, 0.0323059, 0.0167932, 0.0146565, 0.0103066, 0.0115815, 0.00766509, 0.00857287, 0.00540958, 0.00850159,
  1e-99, 1e-99, 1e-99, 0.0407414, 0.028319, 0.0209403, 0.0152207, 0.00921659, 0.00998668, 0.00756254, 0.00715564, 0.00688073,
  1e-99, 1e-99, 1e-99, 0.0318692, 0.0195485, 0.014328, 0.0152939, 0.00760563, 0.00688172, 0.0116874, 0.00670017, 0.0113717,
  1e-99, 1e-99, 1e-99, 0.0458696, 0.0314498, 0.0259211, 0.0173028, 0.0112142, 0.0157959, 0.0122016, 0.01841, 0.0185615,
  1e-99, 1e-99, 1e-99, 0.0405167, 0.0289881, 0.0213844, 0.0131818, 0.0131666, 0.014629, 0.00875657, 0.00840336, 0.0205761,
  1e-99, 1e-99, 1e-99, 0.039329, 0.025446, 0.0176829, 0.0144049, 0.0135103, 0.0148258, 0.0129737, 0.0087925, 0.0186537,
  1e-99, 1e-99, 1e-99, 0.0416085, 0.0253738, 0.0215131, 0.0156886, 0.0119966, 0.0131004, 0.00862069, 0.0153846, 0.0378007,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99
) ),

tauIDFactorizationByPtVsEta_signalAnalysisTauSelectionCaloTauCutBased_CoefficientUncertainty = cms.untracked.vdouble( *(
  1e-99, 1e-99, 1e-99, 0.00222058, 0.00222716, 0.00245066, 0.00262948, 0.00281904, 0.00417362, 0.00352526, 0.00644005, 0.00892248,
  1e-99, 1e-99, 1e-99, 0.00207947, 0.00212347, 0.00223776, 0.00276544, 0.00268412, 0.00305776, 0.00332994, 0.00359496, 0.0102041,
  1e-99, 1e-99, 1e-99, 0.00144792, 0.0014335, 0.00154013, 0.00171341, 0.00211353, 0.00246304, 0.00257013, 0.0026688, 0.00504346,
  1e-99, 1e-99, 1e-99, 0.00162209, 0.00152943, 0.00160119, 0.0018936, 0.00218775, 0.00230374, 0.0020563, 0.00285507, 0.0049554,
  1e-99, 1e-99, 1e-99, 0.00185435, 0.00168361, 0.00192272, 0.00189157, 0.00226443, 0.00320031, 0.00254046, 0.00297482, 0.00411077,
  1e-99, 1e-99, 1e-99, 0.00122245, 0.00123482, 0.00126535, 0.00144938, 0.00155641, 0.00155949, 0.00149869, 0.00211393, 0.00247083,
  1e-99, 1e-99, 1e-99, 0.00161889, 0.00143845, 0.00160934, 0.00183306, 0.00222799, 0.00231481, 0.00230512, 0.0022149, 0.00387335,
  1e-99, 1e-99, 1e-99, 0.00240996, 0.00210669, 0.00240549, 0.00226667, 0.00298954, 0.00400102, 0.00302542, 0.00144509, 0.00527241,
  1e-99, 1e-99, 1e-99, 0.00144676, 0.0013101, 0.001393, 0.00138723, 0.0017845, 0.00220933, 0.0011976, 0.00205094, 0.00201273,
  1e-99, 1e-99, 1e-99, 0.00114371, 0.00100832, 0.00106869, 0.00124321, 0.00119707, 0.00142272, 0.00102973, 0.00132871, 0.0017469,
  1e-99, 1e-99, 1e-99, 0.00158369, 0.00158403, 0.00144393, 0.00178537, 0.00186523, 0.00154321, 0.00161576, 0.0022934, 0.00136426,
  1e-99, 1e-99, 1e-99, 0.00129697, 0.00119719, 0.00124017, 0.00120233, 0.00136771, 0.00144514, 0.00126359, 0.00178599, 0.0017721,
  1e-99, 1e-99, 1e-99, 0.001228, 0.0010975, 0.00113865, 0.00110812, 0.00110963, 0.00155287, 0.0011781, 0.00104881, 0.00163845,
  1e-99, 1e-99, 1e-99, 0.00100544, 0.000906119, 0.000891784, 0.000926893, 0.00116307, 0.00104339, 0.00104933, 0.000905524, 0.00156883,
  1e-99, 1e-99, 1e-99, 0.00119273, 0.0011436, 0.00111233, 0.00114185, 0.00125125, 0.00144718, 0.0013752, 0.00109665, 0.0019685,
  1e-99, 1e-99, 1e-99, 0.00192309, 0.0017071, 0.0018291, 0.00160845, 0.0017539, 0.00180774, 0.00137351, 0.00210084, 0.00377124,
  1e-99, 1e-99, 1e-99, 0.00111085, 0.00105469, 0.00106692, 0.00121043, 0.00131542, 0.00134921, 0.00115312, 0.00124832, 0.00189203,
  1e-99, 1e-99, 1e-99, 0.00229414, 0.00200837, 0.00222316, 0.00249191, 0.0028677, 0.00324051, 0.00204978, 0.00183664, 0.00297603,
  1e-99, 1e-99, 1e-99, 0.0016317, 0.00138981, 0.00160876, 0.00165037, 0.00215062, 0.00212591, 0.00207923, 0.00204463, 0.00300577,
  1e-99, 1e-99, 1e-99, 0.00188731, 0.00191801, 0.00203391, 0.00215253, 0.00206089, 0.00257855, 0.00209747, 0.0025299, 0.00280905,
  1e-99, 1e-99, 1e-99, 0.00123957, 0.00120086, 0.00130255, 0.00167872, 0.0014637, 0.00172043, 0.00206605, 0.00193417, 0.00284293,
  1e-99, 1e-99, 1e-99, 0.00178818, 0.00181576, 0.00205568, 0.00209827, 0.00208243, 0.00309782, 0.00254421, 0.00392503, 0.00464037,
  1e-99, 1e-99, 1e-99, 0.00218135, 0.00229891, 0.00245295, 0.0024478, 0.00294413, 0.00390978, 0.00276907, 0.00343066, 0.00650674,
  1e-99, 1e-99, 1e-99, 0.00121896, 0.00122004, 0.00130716, 0.00147791, 0.00180538, 0.00234416, 0.00210461, 0.00227021, 0.00388956,
  1e-99, 1e-99, 1e-99, 0.00241032, 0.0023976, 0.00277733, 0.00301926, 0.00320622, 0.00436681, 0.00351938, 0.00581484, 0.0113973,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99
) ),

)
