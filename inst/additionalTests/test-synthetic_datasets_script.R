print("SSSTART 1")
c_strDir <- file.path(getwd( ),"..")
print("SSSTART 2")
source(file.path(c_strDir,"synthetic_datasets_script.R"))

########## General Controls
# Control the percision used for rounding in the test results.
c_iRoundingPercision = 2
c_fVerbose = TRUE
sCalibrateFile = "Calibrate-RLNorm.tsv"
sCalibrateSparseFile = "Calibrate-Sparse.tsv"

########## Original calibration file data
c_vdExpLog = c( -0.877070018720874,1.24126858906963,-0.71743987312899,-0.579818495252942,-0.916290731874155,-1.42711635564015,-0.350976922824095,-2.68824757380603,4.45611317802658,-1.39432653281715,-0.601479992034121,-1.03282454813011,4.6737256284836,2.52348613175548,-0.867500567704723,-1.6094379124341,-2.57702193869581,0.799307376388336,1.99768903980758,-1.58963528513792,7.55066334615255,3.13948623719869,1.50318811259139,0.0353671438372913,2.43220886078755,4.85881966369566,1.67372640231646,1.92629036218566,3.4618533755355,1.93036131866568,4.68792920115586,-0.0921152889078056,1.84498423046535,-2.02495335639577,-1.83258146374831,2.14710019015365,-0.785262469467751,1.0123279200711,0.542324290825362,3.30630044979247,-1.65025990695436,0.525911261184032,2.43571640559723,-0.653926467406664,4.43585172136879,0.790273891290668,0.72270598280149,0.0506931143155182,-0.71743987312899,1.21669157673371,5.10315370174078,1.5295285292058,-1.20397280432594,3.88362353090645,-2.52572864430826,1.82841278687609,-1.20397280432594,4.13523055469444,-0.524248644098131,-0.551647618286246,-0.0790432073404529,-1.24479479884619,6.32837954283977,-1.10262031006565,-0.20334092401803,-0.415515443961666,0.5607579925142,-0.233193887167711,-2.52572864430826,5.19168938132325,5.3469170363841,5.12508240135524,-1.80788885115794,-1.65025990695436,5.59600422381589,6.02419307773006,0.845009529869191,-2.52572864430826,-0.317454230785451,2.86539577182599,2.59614982860398,-0.776528789498996,0.501986675098786,-0.742337424750717,0.249980205267769,3.63989927189222,-1.09064411901893,6.50169199200066,2.4308023907948,2.11529119457353,4.66317490827467,1.17803942229943,4.3994244118503,5.69332801675234,3.43501893013767,-3.03655426807425,-0.988861424708991,4.46558622777281,2.94528073005789,3.11688739511046,4.15211019907224,3.28929700916393,2.14006616349627,5.78226088252746,0.681074599325676,3.17588481306311,5.01359092192189,1.95103982687531,0.877134016672961,2.47754629538171,4.2285840374303,4.32772883126603,-1.57021719928082,-2.88240358824699,-1.96611285637283,-1.78379129957888,0.29266961396282,2.67662833109387,1.04802050255205,1.04942204447734,-1.58963528513792,0.0544881852840698,-0.0877389143080067,1.38729386145297,-2.52572864430826,-0.2484613592985,-1.80788885115794,-0.867500567704723,-3.12356564506388,-1.76026080216868,-0.988861424708991,-0.877070018720874,0.289680075114454,-0.53785429615391,-2.15416508787577,-1.99510039324608,0.666803205220343,1.70765295993106,-2.52572864430826,5.39355481643488,1.92861865194525,1.27759494419655,0.227932068046007,-2.12026353620009,4.45066611039057,0.0769610411361284,1.26186428274171,1.57608793275255,3.5184469417124,0.710987098688276,-2.95651156040071,0.4265740713184,1.86128553187667,-0.661648513500574,-2.12026353620009,2.77832225408754,5.86404018285232,5.91260530552877 )

c_vdExp = c(0.416,3.46,0.488,0.56,0.4,0.24,0.704,0.068,86.152,0.248,0.548,0.356,107.096,12.472,0.42,0.2,0.076,2.224,7.372,0.204,1902.004,23.092,4.496,1.036,11.384,128.872,5.332,6.864,31.876,6.892,108.628,0.912,6.328,0.132,0.16,8.56,0.456,2.752,1.72,27.284,0.192,1.692,11.424,0.52,84.424,2.204,2.06,1.052,0.488,3.376,164.54,4.616,0.3,48.6,0.08,6.224,0.3,62.504,0.592,0.576,0.924,0.288,560.248,0.332,0.816,0.66,1.752,0.792,0.08,179.772,209.96,168.188,0.164,0.192,269.348,413.308,2.328,0.08,0.728,17.556,13.412,0.46,1.652,0.476,1.284,38.088,0.336,666.268,11.368,8.292,105.972,3.248,81.404,296.88,31.032,0.048,0.372,86.972,19.016,22.576,63.568,26.824,8.5,324.492,1.976,23.948,150.444,7.036,2.404,11.912,68.62,75.772,0.208,0.056,0.14,0.168,1.34,14.536,2.852,2.856,0.204,1.056,0.916,4.004,0.08,0.78,0.164,0.42,0.044,0.172,0.372,0.416,1.336,0.584,0.116,0.136,1.948,5.516,0.08,219.984,6.88,3.588,1.256,0.12,85.684,1.08,3.532,4.836,33.732,2.036,0.052,1.532,6.432,0.516,0.12,16.092,352.144,369.668)

c_vdPercentZero = c( 0.964,0.672,0.988,0.92,0.988,0.952,0.88,0.984,0.428,0.936,0.864,0.916,0.472,0.508,0.972,0.956,0.98,0.92,0.672,0.952,0.064,0.72,0.832,0.984,0.66,0.332,0.864,0.908,0.728,0.884,0.392,0.888,0.748,0.96,0.976,0.84,0.952,0.972,0.788,0.948,0.956,0.784,0.848,0.944,0.7,0.972,0.9,0.884,0.924,0.776,0.096,0.696,0.936,0.58,0.984,0.772,0.936,0.352,0.96,0.88,0.932,0.936,0.072,0.952,0.92,0.912,0.876,0.92,0.968,0.136,0.052,0.12,0.988,0.944,0.2,0.02,0.88,0.98,0.916,0.464,0.452,0.892,0.8,0.932,0.78,0.22,0.936,0.212,0.504,0.608,0.484,0.756,0.256,0.064,0.752,0.988,0.944,0.54,0.948,0.824,0.556,0.492,0.548,0.016,0.788,0.932,0.092,0.48,0.784,0.396,0.152,0.732,0.972,0.984,0.98,0.98,0.932,0.7,0.884,0.972,0.924,0.916,0.936,0.94,0.984,0.944,0.964,0.936,0.984,0.988,0.964,0.98,0.912,0.952,0.972,0.972,0.928,0.728,0.976,0.416,0.932,0.952,0.952,0.984,0.416,0.876,0.968,0.772,0.608,0.82,0.988,0.948,0.892,0.932,0.98,0.84,0.216,0.612 )

c_vdSD = c(1.29288222850903,1.11656790017514,1.92641197859569,1.02105677634939,1.48524719110303,0.870378434477115,0.959635078912515,0.583945419522526,1.79348082582133,0.669655926119326,0.600623613084851,0.853648651939737,1.72677780263148,1.26147280133706,0.734982664869215,0.914826396438766,0.719267497333394,1.36445803823611,1.24480603198208,0.554831763342855,1.70747836197823,1.49834967389745,1.10454940596015,1.86473285616939,1.1385848782746,1.4834337653741,1.52109111483939,1.46795133122398,1.98185399587786,1.45950174080202,1.59366893620959,1.03774612854838,1.44800233026963,0.592107847055163,0.594358593206619,1.6320293719476,0.962434020985145,2.07513585839894,1.09851647985654,2.49835214683079,0.976161247538305,1.10000715772956,1.80999152731209,1.15455694103497,2.29447750308131,1.84856812708444,1.18034213829724,0.938036155588448,0.874004621310478,1.09917849384688,1.7434683468884,1.16454079297455,0.728859758118786,1.77099312319428,0.694335857646383,1.37928702782278,0.76009658704113,1.68201745534149,1.35054907979964,0.832670104071402,1.16865528891661,0.847430958060489,1.5524935357744,0.904423733977812,1.17860611568599,1.03253595130507,1.29763235709243,1.06536965945072,0.27424593232245,1.57376990950369,1.34803070031666,1.55956753125577,1.10470442660092,0.723280490439043,1.5977407349549,0.930850462498865,1.29438882062999,0.803822732571653,0.986243910474401,1.29553952694038,1.28823184910116,0.685283831641858,0.946823554040535,0.918659788437422,0.788068631230099,1.22250675713921,0.819552326795716,1.59822401385881,1.17461364613988,1.17170376599601,1.87008425727893,1.11396621015164,1.41133141289528,1.38435992039007,1.37755755926622,0.256091402042052,0.758162389852548,1.64572542671349,2.14130589750552,1.88179351330357,1.41089291208342,1.60224644297979,1.14667348092866,0.942793883548745,0.922341311891637,1.96208871494235,1.56299028564406,1.01027693685149,1.02567334310102,1.17471546316181,1.46981919823393,2.08186319151382,1.17306243431671,0.758804728095711,0.843696416029418,0.670981308304214,1.15201957724894,1.56952695767858,1.42621506263831,2.07361684190959,0.751160258829315,1.2533916605354,1.41373507230623,1.73889522835292,1.00489987867617,1.23624923115956,0.556618232166134,1.1257008349353,0.612456251969496,0.824857904639158,1.25061326841638,1.27907129252095,1.03099622126377,1.1743855423561,0.685788698475515,1.22474258547412,1.49553187222792,1.23068112783264,0.421237925499913,2.25088768513035,2.24179750266376,1.88939818126567,1.42078414973316,0.830682791171414,1.8296524259332,1.16414627194934,2.22431472689055,1.27214629360624,1.64269420055488,1.20601967125448,0.294925311390253,1.7163875280654,1.62333836869266,1.16870601697591,0.618214596868589,1.751366359945,2.40861367145808,1.4507511261923)

########## Start testing


### Test that the script runs to completion
#!#

### funcBoxPlotOutliers
### Plotting function.


### funcCalibrateRLNormToMicrobiome
### TODO(Timothy): Make a file and calculate the correct relationships and make sure that is found here.


### funcEstimateFeatureSD
### Pass 10-22-2013
context("Test funcEstimateFeatureSD")
test_that( "Test that funcEstimateFeatureSD gives the correct answer.", {
  expect_equal( 0, funcEstimateFeatureSD( 0, 12834 ) )
  expect_equal( 1, funcEstimateFeatureSD( 1, 1 ) )
  expect_equal( 1.5, funcEstimateFeatureSD( .5, 3 ) )
  expect_equal( 1010, funcEstimateFeatureSD( 101, 10 ) )
  expect_equal( -24, funcEstimateFeatureSD( -3, 8 ) )
})


### funcEstimateGrandSD
### PASS 10-22-2013
context("Test funcEstimateGrandSD")
test_that( "Test that funcEstimateGrandSD gives the correct answer.", {
  expect_equal( 0, funcEstimateGrandSD( 0, 12834 ) )
  expect_equal( 1, funcEstimateGrandSD( 1, 1 ) )
  expect_equal( 1.5, funcEstimateGrandSD( .5, 3 ) )
  expect_equal( 1010, funcEstimateGrandSD( 101, 10 ) )
  expect_equal( -24, funcEstimateGrandSD( -3, 8 ) )
})


### funcEstimatePercentZero
### PASS 10-22-2013
context("Test funcEstimatePercentZero")
test_that( "Test that funcEstimatePercentZero gives the correct answer.", {
  expect_equal( 0.125, funcEstimatePercentZero( vdExpLog = .5, dBetaZero = .2, dBeta2Zero = .1 ) )
  expect_equal( 0.125, funcEstimatePercentZero( vdExpLog = .5, dBetaZero = .2, dBeta2Zero = .1, dInterceptZero = 0 ) )
  expect_equal( 0.125, funcEstimatePercentZero( vdExpLog = .5, dBetaZero = .2, dBeta2Zero = .1, dInterceptZero = 0, dScale = 1 ) )
  expect_equal( 0.39, funcEstimatePercentZero( vdExpLog = .3, dBetaZero = 1, dBeta2Zero = 1, dInterceptZero = 0, dScale = 1 ) )
  expect_equal( 1, funcEstimatePercentZero( vdExpLog = .3, dBetaZero = 5, dBeta2Zero = 1.3, dInterceptZero = 3, dScale = 2 ) )
  expect_equal( 0, funcEstimatePercentZero( vdExpLog = .67, dBetaZero = .3, dBeta2Zero = .1, dInterceptZero = -77, dScale = .2 ) )
})


###funcForceMinCountsInMinSamples
# PASS 10-22-2013
vdFeature1 = rep(1,10)
vdFeature2 = 1:10
context("Test funcForceMinCountsInMinSamples")
test_that( "Test that funcForceMinCountsInMinSamples increases features to the min.",
{
  expect_equal( 10, length( which( funcForceMinCountsInMinSamples(vdFeature1, 5, 10) >= 5 ) ) )
  expect_equal( 5, length( which( funcForceMinCountsInMinSamples(vdFeature1, 5, 5) >= 5 ) ) )
  expect_equal( 10, length( which( funcForceMinCountsInMinSamples(vdFeature1, 20, 10) >= 20 ) ) )
  expect_equal( 10, length( which( funcForceMinCountsInMinSamples(vdFeature2, 5, 10) >= 5 ) ) )
  expect_equal( 6, length( which( funcForceMinCountsInMinSamples(vdFeature2, 5, 5) >= 5 ) ) )
  expect_equal( 10, length( which( funcForceMinCountsInMinSamples(vdFeature2, 20, 10) >= 20 ) ) )
  expect_equal( vdFeature1, funcForceMinCountsInMinSamples(vdFeature1, 0, 10) )
  expect_equal( vdFeature2, funcForceMinCountsInMinSamples(vdFeature2, 1, 10) )
  expect_equal( vdFeature2, funcForceMinCountsInMinSamples(vdFeature2, 0, 0) )
})


### funcGenerateExpVectorParameters
### PASS 10-22-2013
context( "Test funcGenerateExpVectorParameters" )
dMu <- 2
dSD <- 3
iCount = 1000000
dLogMu <- log( dMu )
dLogSD <- log( dSD )
dExp <- funcGetExp( dMu, dSD )
dRelationship <- dLogSD / dExp
dReadDepth <- dExp * iCount
vdDistribution <- rlnorm( iCount, dLogMu, dLogSD )
lResult = funcGenerateExpVectorParameters(vdDistribution)
test_that( "Test that funcGenerateExpVectorParameters gives the correct answer for good case 1.", {
  expect_equal( round( dLogMu, 1 ), round( lResult$GrandLogMu, 1 ) )
  expect_equal( round( dExp, 1), round( lResult$GrandExp, 1 ) )
  expect_equal( round( dLogSD, 1), round( lResult$GrandLogSD, 1 ) )
  expect_true( dReadDepth - ( dReadDepth * 0.1 ) < lResult$dReadDepth )
  expect_true( dReadDepth + ( dReadDepth * 0.1 ) > lResult$dReadDepth )
  expect_equal( round( dRelationship, 1 ), round( lResult$GrandLogSDBeta, 1 ) )
})

dMu <- 5
dSD <- 4
iCount = 1000000
dLogMu <- log( dMu )
dLogSD <- log( dSD )
dExp <- funcGetExp( dMu, dSD )
dRelationship <- dLogSD / dExp
dReadDepth <- dExp * iCount
vdDistribution <- rlnorm( iCount, dLogMu, dLogSD )
lResult = funcGenerateExpVectorParameters(vdDistribution)
test_that( "Test that funcGenerateExpVectorParameters gives the correct answer for good case 2.", {
  expect_equal( round( dLogMu, 1 ), round( lResult$GrandLogMu, 1 ) )
  expect_equal( round( dExp, 1), round( lResult$GrandExp, 1 ) )
  expect_equal( round( dLogSD, 1), round( lResult$GrandLogSD, 1 ) )
  expect_true( dReadDepth - ( dReadDepth * 0.1 ) < lResult$dReadDepth )
  expect_true( dReadDepth + ( dReadDepth * 0.1 ) > lResult$dReadDepth )
  expect_equal( round( dRelationship, 1 ), round( lResult$GrandLogSDBeta, 1 ) )
})


### func_generate_bug_bug_spiking_matrix
#!#


### func_generate_lefse_matrices
#!#


### funcGenerateFeatureParameters
context( "Test funcGenerateFeatureParameters" )
iNumberFeatures <- 10
iNumberSamples <- 1
iMinNumberSamples <- 0
iReadDepth <- 365
vdExp <- NA
vdMu <- NA
vdSD <- NA
vdPercentZero <- NA
dBetaGrandSD <- log( 3 ) / ( iReadDepth / iNumberFeatures )
lSDRel = list(BetaSD = 1, InterceptSD = 2 )
lPercentZeroRel = list(InterceptZero = 2.3, BetaZero = 0, Beta2Zero = 1.56)

lResult = funcGenerateFeatureParameters( int_number_features = iNumberFeatures, int_number_samples = iNumberSamples, iMinNumberSamples = iMinNumberSamples, iReadDepth = iReadDepth, vdExp = vdExp, vdMu = vdMu, vdSD = vdSD, vdPercentZero = vdPercentZero, lSDRel = lSDRel, lPercentZeroRel = lPercentZeroRel, dBetaGrandSD = dBetaGrandSD )
test_that( "funcGenerateFeatureParameters: Generating all vectors.", {
  expect_equal( round( iReadDepth / iNumberFeatures ), round( mean( lResult$exp ) ) )
  expect_equal( iReadDepth, ceiling( sum( sapply( 1:length( lResult$exp ), function( x ) funcGetExp( lResult$mu[ x ], lResult$sd[ x ]) ) ) ) )
  expect_equal( lResult$exp, sapply( 1:length( lResult$mu ), function(x) funcGetExp( lResult$mu[x], lResult$sd[x] ) ) )
  expect_true( length( lResult$PercentZero ) == 0 )
})

iNumberFeatures <- 10
iNumberSamples <- 1
iMinNumberSamples <- 0
iReadDepth <- 365
vdExp <- c(37.45276, 40.38829, 29.10671, 41.38829, 38.95904, 33.30842, 37.38829, 30.38829, 34.38829, 41.38829)
vdMu <- NA
vdSD <- NA
vdPercentZero <- NA
dBetaGrandSD <- log( 3 ) / ( iReadDepth / iNumberFeatures )
lSDRel = list(BetaSD = 1.1, InterceptSD = 7.2 )
lPercentZeroRel = list(InterceptZero = .32, BetaZero = 2, Beta2Zero = .56)

lResult = funcGenerateFeatureParameters( int_number_features = iNumberFeatures, int_number_samples = iNumberSamples, iMinNumberSamples = iMinNumberSamples, iReadDepth = iReadDepth, vdExp = vdExp, vdMu = vdMu, vdSD = vdSD, vdPercentZero = vdPercentZero, lSDRel = lSDRel, lPercentZeroRel = lPercentZeroRel, dBetaGrandSD = dBetaGrandSD )
test_that( "funcGenerateFeatureParameters: Starting with Mu and SD vector.", {
  expect_equal( round( iReadDepth / iNumberFeatures ), round( mean( lResult$exp ) ) )
  expect_equal( iReadDepth, ceiling( sum( sapply( 1:length( lResult$exp ), function( x ) funcGetExp( lResult$mu[ x ], lResult$sd[ x ]) ) ) ) )
  expect_equal( lResult$exp, sapply( 1:length( lResult$mu ), function(x) funcGetExp( lResult$mu[x], lResult$sd[x] ) ) )
#TODO#
  expect_equal( 0, sum( lResult$PercentZero ) )
  expect_equal( round( lResult$mu, 2), round( c(0.05285117, 0.04323700, 0.09918717, 0.04046334, 0.04762098, 0.07139393, 0.05309049, 0.08947141, 0.06587467, 0.04046334), 2 ) )
  expect_equal( round( lResult$sd, 2), round( c(37.45276, 40.38829, 29.10671, 41.38829, 38.95904, 33.30842, 37.38829, 30.38829, 34.38829, 41.38829), 2 ) )
  expect_equal( round( lResult$PercentZero, 2), round( rep( 1, length( lResult$exp ) ), 2 ) )
})

iNumberFeatures <- 10
iNumberSamples <- 1
iMinNumberSamples <- 0
iReadDepth <- 365
vdExp <- c(37.45276, 40.38829, 29.10671, 41.38829, 38.95904, 33.30842, 37.38829, 30.38829, 34.38829, 41.38829)
vdMu <- c(0.05285117, 0.04323700, 0.09918717, 0.04046334, 0.04762098, 0.07139393, 0.05309049, 0.08947141, 0.06587467, 0.04046334)
vdSD <- NA
vdPercentZero <- NA
lSDRel = list(BetaSD = 5, InterceptSD = 8.4 )
lPercentZeroRel = list(InterceptZero = 2.1, BetaZero = 41.2, Beta2Zero = 1.756)
dBetaGrandSD <- log( 3 ) / ( iReadDepth / iNumberFeatures )

lResult = funcGenerateFeatureParameters( int_number_features = iNumberFeatures, int_number_samples = iNumberSamples, iMinNumberSamples = iMinNumberSamples, iReadDepth = iReadDepth, vdExp = vdExp, vdMu = vdMu, vdSD = vdSD, vdPercentZero = vdPercentZero, lSDRel = lSDRel, lPercentZeroRel = lPercentZeroRel, dBetaGrandSD = dBetaGrandSD )
test_that( "funcGenerateFeatureParameters: Starting with Mu and SD and Mu vector.", {
  expect_equal( round( iReadDepth / iNumberFeatures, 2 ), round( mean( lResult$exp ), 2 ) )
  expect_equal( iReadDepth, sum( sapply( 1:length( lResult$exp ), function( x ) funcGetExp( lResult$mu[ x ], lResult$sd[ x ]) ) ) )
  expect_equal( round( lResult$exp, 2), round( lResult$sd, 2) )
  expect_equal( 0, sum( lResult$PercentZero ) )
  expect_equal( round( lResult$sd, 2), round( c(37.45276, 40.38829, 29.10671, 41.38829, 38.95904, 33.30842, 37.38829, 30.38829, 34.38829, 41.38829), 2 ) )
  expect_equal( round( lResult$PercentZero, 2), round( rep( 1, length( lResult$exp ) ), 2 ) )
})

iNumberFeatures <- 10
iNumberSamples <- 1
iMinNumberSamples <- 0
iReadDepth <- 365
vdExp <- c(37.45276, 40.38829, 29.10671, 41.38829, 38.95904, 33.30842, 37.38829, 30.38829, 34.38829, 41.38829)
vdMu <- c(0.05285117, 0.04323700, 0.09918717, 0.04046334, 0.04762098, 0.07139393, 0.05309049, 0.08947141, 0.06587467, 0.04046334)
vdSD <- c(37.45276, 40.38829, 29.10671, 41.38829, 38.95904, 33.30842, 37.38829, 30.38829, 34.38829, 41.38829)
vdPercentZero <- NA
lSDRel = list(BetaSD = 1, InterceptSD = 2 )
lPercentZeroRel = list(InterceptZero = 2.3, BetaZero = 0, Beta2Zero = 1.56)
dBetaGrandSD <- log( 3 ) / ( iReadDepth / iNumberFeatures )

lResult = funcGenerateFeatureParameters( int_number_features = iNumberFeatures, int_number_samples = iNumberSamples, iMinNumberSamples = iMinNumberSamples, iReadDepth = iReadDepth, vdExp = vdExp, vdMu = vdMu, vdSD = vdSD, vdPercentZero = vdPercentZero, lSDRel = lSDRel, lPercentZeroRel = lPercentZeroRel, dBetaGrandSD = dBetaGrandSD )
test_that( "funcGenerateFeatureParameters: Starting with Mu and SD and Mu vector.", {
  expect_equal( round( iReadDepth / iNumberFeatures, 2 ), round( mean( lResult$exp ), 2 ) )
  expect_equal( iReadDepth, sum( sapply( 1:length( lResult$exp ), function( x ) funcGetExp( lResult$mu[ x ], lResult$sd[ x ]) ) ) )
  expect_equal( round( lResult$exp, 2), round( lResult$sd, 2) )
  expect_equal( 0, sum( lResult$PercentZero ) )
  expect_equal( round( lResult$sd, 2), round( c(37.45276, 40.38829, 29.10671, 41.38829, 38.95904, 33.30842, 37.38829, 30.38829, 34.38829, 41.38829), 2 ) )
  expect_equal( round( lResult$PercentZero, 2), round( rep( 1, length( lResult$exp ) ), 2 ) )
})

iNumberFeatures <- 10
iNumberSamples <- 1
iMinNumberSamples <- 0
iReadDepth <- 365
vdExp <- c(37.45276, 40.38829, 29.10671, 41.38829, 38.95904, 33.30842, 37.38829, 30.38829, 34.38829, 41.38829)
vdMu <- c(0.05285117, 0.04323700, 0.09918717, 0.04046334, 0.04762098, 0.07139393, 0.05309049, 0.08947141, 0.06587467, 0.04046334)
vdSD <- c(37.45276, 40.38829, 29.10671, 41.38829, 38.95904, 33.30842, 37.38829, 30.38829, 34.38829, 41.38829)
vdPercentZero <- rep( 1, length( lResult$exp ) )
lSDRel = list(BetaSD = 1, InterceptSD = 2 )
lPercentZeroRel = list(InterceptZero = 2.3, BetaZero = 1, Beta2Zero = 1.56)
dBetaGrandSD <- log( 3 ) / ( iReadDepth / iNumberFeatures )

lResult = funcGenerateFeatureParameters( int_number_features = iNumberFeatures, int_number_samples = iNumberSamples, iMinNumberSamples = iMinNumberSamples, iReadDepth = iReadDepth, vdExp = vdExp, vdMu = vdMu, vdSD = vdSD, vdPercentZero = vdPercentZero, lSDRel = lSDRel, lPercentZeroRel = lPercentZeroRel, dBetaGrandSD = dBetaGrandSD )
test_that( "funcGenerateFeatureParameters: Starting with Mu and SD and Mu and PZ vector.", {
  expect_equal( round( iReadDepth / iNumberFeatures, 2 ), round( mean( lResult$exp ), 2 ) )
  expect_equal( iReadDepth, sum( sapply( 1:length( lResult$exp ), function( x ) funcGetExp( lResult$mu[ x ], lResult$sd[ x ]) ) ) )
  expect_equal( round( lResult$exp, 2), round( lResult$sd, 2) )
  expect_equal( 0, sum( lResult$PercentZero ) )
})

iNumberFeatures <- 100
lResult = funcGenerateFeatureParameters( int_number_features = iNumberFeatures, int_number_samples = iNumberSamples, iMinNumberSamples = iMinNumberSamples, iReadDepth = iReadDepth, vdExp = vdExp, vdMu = vdMu, vdSD = vdSD, vdPercentZero = vdPercentZero, lSDRel = lSDRel, lPercentZeroRel = lPercentZeroRel, dBetaGrandSD = dBetaGrandSD )
testthat( "funcGenerateFeatureParameters: Starting with all vectors but wrong size.", {
  expect_equal(,)
})

test_that( "funcGenerateFeatureParameters: Starting with min number really high, min samples and sd = 1.", {
  expect_equal(iNumberFeatures,length(lResult$mu))
  expect_equal(iNumberFeatures,length(lResult$exp))
  expect_equal(iNumberFeatures,length(lResult$sd))
  expect_equal(iNumberFeatures,length(lResult$PercentZero))
})

# iNumberFeatures <- 1000
# iNumberSamples <- 1
# iMinNumberSamples <- 0
# iReadDepth <- 36500
# vdExp <- NA
# vdMu <- NA
# vdSD <- NA
# vdPercentZero <- NA
# dBetaSD <- 1
# dBetaZero <- 0
# dBetaGrandSD <- log( 3 ) / ( iReadDepth / iNumberFeatures )
# 
# testthat( "funcGenerateFeatureParameters: Read depth 36500." {
#   expect_equal(36500,sum(lResult$exp))
# })


### func_generate_metadata
#!#


### func_generate_random_lognormal_matrix
#!#


### func_generate_random_lognormal_with_outliers
context("Test func_generate_random_lognormal_with_outliers")
#!#


### func_generate_random_lognormal_with_multivariate_spikes
context("Test func_generate_random_lognormal_with_multivariate_spikes")
#!#


### funcGetExp
# PASS 10-22-2013
context("Test funcGetExp")
dSDOne = round(3,4)
dMuOne = round(4,4)
dExpOne = mean(rlnorm(1000000,log(dMuOne),log(dSDOne)))
dSDTwo = round(2,4)
dMuTwo = round(3,4)
dExpTwo = mean(rlnorm(1000000,log(dMuTwo),log(dSDTwo)))
dSDThree = round(1,4)
dMuThree = round(1,4)
dExpThree = mean(rlnorm(1000000,log(dMuThree),log(dSDThree)))
test_that("Test that the expected SD is calculated correctly.",{
  expect_equal(round(dExpOne,1),round(funcGetExp(dMuOne,dSDOne),1))
  expect_equal(round(dExpTwo,1),round(funcGetExp(dMuTwo,dSDTwo),1))
  expect_equal(round(dExpThree,1),round(funcGetExp(dMuThree,dSDThree),1))
})


### funcGetMu
# PASS 10-22-2013
context("Test funcGetMu")
dSDOne = round(3,4)
dMuOne = round(4,4)
dExpOne = mean(rlnorm(1000000,log(dMuOne),log(dSDOne)))
dSDTwo = round(2,4)
dMuTwo = round(3,4)
dExpTwo = mean(rlnorm(1000000,log(dMuTwo),log(dSDTwo)))
dSDThree = round(1,4)
dMuThree = round(1,4)
dExpThree = mean(rlnorm(1000000,log(dMuThree),log(dSDThree)))
test_that("Test that the expected SD is calculated correctly.",{
  expect_equal(dMuOne,round(funcGetMu(dExpOne,dSDOne),1))
  expect_equal(dMuTwo,round(funcGetMu(dExpTwo,dSDTwo),1))
  expect_equal(dMuThree,round(funcGetMu(dExpThree,dSDThree),1))
})


### funcGetParamsForReadDepth
# PASS 10-22-2013
context("Test funcGetParamsForReadDepth")
dMu <- 2
dSD <- 3
dExp <- funcGetExp( dMu, dSD )
iCount <- 1000000
vdDistribution <- rlnorm( iCount, log( dMu ), log( dSD ) )
lResult <- funcGetParamsForReadDepth( mean( vdDistribution ), log(dSD) / dExp )

test_that( "Test that the parameters measured by funcGetParamsForExpectation are off less than 1 decimal, test 1.", {
  expect_equal( round( log(dMu) ), round( lResult$dLogMu ) )
  expect_equal( round( log(dSD) ), round( lResult$dLogSD ) )
})

dMu <- 4
dSD <- 2
dExp <- funcGetExp( dMu, dSD )
iCount <- 1000000
vdDistribution <- rlnorm( iCount, log( dMu ), log( dSD ) )
lResult <- funcGetParamsForReadDepth( mean( vdDistribution ), log(dSD) / dExp )

test_that( "Test that the parameters measured by funcGetParamsForExpectation are off less than 1 decimal, test 1.", {
  expect_equal( round( log(dMu) ), round( lResult$dLogMu ) )
  expect_equal( round( log(dSD) ), round( lResult$dLogSD ) )
})


### funcGetRowMetric
# PASS 10-22-2013
viOne <- 1:10
viTwo <- c( 1, 5, 2, 8, 3, 4, 6, 9, 12, 2 )
viThree <- c( 1, 24, 32, 5, 21, 2, 4, 6, 8, 67 )
vdAnswer3 <- c( 5.5, 5.2, 17 )
test_that( "Test that the mean is calculated correctly for vectors or matrices.", {
  expect_equal( funcGetRowMetric( viOne, mean ), 5.5 )
  expect_equal( funcGetRowMetric( matrix( viOne, ncol = length( viOne ), byrow = TRUE ), mean ), 5.5 )
  expect_equal( funcGetRowMetric( matrix( c( viOne, viTwo, viThree ), ncol = length( viOne ), byrow = TRUE ), mean), vdAnswer3 )
  expect_equal( funcGetRowMetric( viOne, sum ), 55 )
  expect_equal( funcGetRowMetric( matrix( viOne, ncol = length( viOne ), byrow = TRUE ), sum ), 55 )
})


### Test funcGetSD
# PASS 10-22-2013
context("Test funcGetSD")
viOne = 1:10
viTwo = c(1,5,2,8,3,4,6,9,12,2)
viThree = c(1,24,32,5,21,2,4,6,8,67)
dSDOne = 3.02765
dSDTwo = 3.552777
dSDThree = 20.51016
test_that("Test that the SD is calculated correctly for vectors or matrices.",{
  expect_equal(round(funcGetSD(viOne),c_iRoundingPercision), round(dSDOne,c_iRoundingPercision))
  expect_equal(round(funcGetSD(matrix(c(viOne),ncol=length(viOne),byrow=TRUE)),c_iRoundingPercision), round(dSDOne,c_iRoundingPercision))
  expect_equal(round(funcGetSD(matrix(c(viOne, viTwo, viThree),ncol=length(viOne),byrow=TRUE)),c_iRoundingPercision), round(c(dSDOne ,dSDTwo, dSDThree),c_iRoundingPercision))
})


### funcIsFactorMetadataValid
# PASS 10-22-2013
context("Test funcIsFactorMetadataValid")
vxMetadataOne = c(1,2,2,1,1,1,2,2,1,2,1,2)
vxMetadataTwo = c(1,2,3,1,1,1,3,3,2,2,3,2,1,3)
test_that("Test that minimal levels of data are detected correctly, discontinuous data only.",{
  expect_equal(funcIsFactorMetadataValid(vxMetadataOne,1), TRUE)
  expect_equal(funcIsFactorMetadataValid(vxMetadataOne,8), FALSE)
  expect_equal(funcIsFactorMetadataValid(vxMetadataTwo,4), TRUE)
  expect_equal(funcIsFactorMetadataValid(vxMetadataTwo,5), FALSE)
})


# Could test more
### funcMakeFeature
# PASS 10-22-2013
context("Test funcMakeFeature for min counts")
# Test for minimum count in a sample
dMuOne = 1
dSDOne = 1
dPercentZeroOne = 0.0
iNumberSamplesOne = 100
iMinNumberCountsOne = 10
iMinNumberSamplesOne = 10
dTotalReadDepth = 1

iMinNumberCountsTwo = 100
iMinNumberCountsThree = 0

test_that("Test that the expected number of samples are given the minimum value.",{
  expect_equal(length(which(funcMakeFeature(dMuOne,dSDOne,dPercentZeroOne,iNumberSamplesOne,iMinNumberCountsOne,iMinNumberSamplesOne, dTotalReadDepth)$Feature>=iMinNumberCountsOne)),iMinNumberSamplesOne)
  expect_equal(length(which(funcMakeFeature(dMuOne,dSDOne,dPercentZeroOne,iNumberSamplesOne,iMinNumberCountsTwo,iMinNumberSamplesOne, dTotalReadDepth)$Feature>=iMinNumberCountsTwo)),iMinNumberSamplesOne)
  expect_equal(length(which(funcMakeFeature(dMuOne,dSDOne,dPercentZeroOne,iNumberSamplesOne,iMinNumberCountsThree,iMinNumberSamplesOne, dTotalReadDepth)$Feature>=iMinNumberCountsThree)), iNumberSamplesOne)
})


### funcNormalizeMicrobiome
# PASS 10-22-2013
context("Test funcNormalizeMicrobiome")
mtrxCountsOne = matrix(c(1,1,1,2,2,2,3,3,3,4,4,4), byrow=TRUE, nrow=4, ncol=3)
mtrxRelAbOne = matrix(c(0.1,0.1,0.1,0.2,0.2,0.2,0.3,0.3,0.3,0.4,0.4,0.4), byrow=TRUE, nrow=4, ncol=3)
mtrxCountsTwo = matrix(c(1,1,1,2,2,2,3,3,3,4,4,4), nrow=4, ncol=3)
mtrxRelAbTwo = matrix(c(1/5,1/5,1/5,2/5,2/10,2/10,3/10,3/10,3/15,4/15,4/15,4/15), nrow=4, ncol=3)
mtrxCountsZero = matrix(rep(0,4*3), byrow=TRUE, nrow=4, ncol=3)
mtrxRelAbZero = matrix(rep(0,4*3), byrow=TRUE, nrow=4, ncol=3)
mtrxCountsFour = matrix(c(1,0,1,0,2,0,3,0,3,0,4,0), byrow=TRUE, nrow=4, ncol=3)
mtrxRelAbFour = matrix(c(1/4,0/6,1/4,0/4,2/6,0/4,3/4,0/6,3/4,0/4,4/6,0/4), byrow=TRUE, nrow=4, ncol=3)
test_that("Test that normalization occurs by column and when zeros are present.",{
  expect_equal(mtrxRelAbOne,funcNormalizeMicrobiome(mtrxCountsOne))
  expect_equal(mtrxRelAbTwo,funcNormalizeMicrobiome(mtrxCountsTwo))
  expect_equal(mtrxRelAbZero,funcNormalizeMicrobiome(mtrxCountsZero))
  expect_equal(mtrxRelAbFour,funcNormalizeMicrobiome(mtrxCountsFour))
})


### funcNumericIsInt
# PASS 10-22-2013
context("Test funcNumericIsInt")
test_that("Test that integers are correctly identified.",{
  expect_equal(funcNumericIsInt(1), TRUE)
  expect_equal(funcNumericIsInt(2.0), TRUE)
  expect_equal(funcNumericIsInt(1.5), FALSE)
  expect_equal(funcNumericIsInt(0.00005), FALSE)
})


### funcPlotSpikeins
### Plotting function


### funcPrepareMetadata
# PASS 10-22-2013
context("Test funcPrepareMetadata with dummying")
vstrMetadataOne = c(1,2,2,1,1,1,2,2,1,2,1,2)
vstrMetadataTwo = c(1,2,3,1,1,3,2,3,1,2,1,2)
vstrMetadataThree = c(1,2,3,4,1,3,2,4,1,2,1,2)

lResult1 = funcPrepareMetadata(c(3),vstrMetadataOne,TRUE)
lResult2 = funcPrepareMetadata(c(8),vstrMetadataOne,TRUE)
lResult3 = funcPrepareMetadata(c(3,13,23),matrix(c(vstrMetadataOne, vstrMetadataTwo, vstrMetadataThree),ncol=length(vstrMetadataThree),byrow=TRUE),TRUE)
lResult4 = funcPrepareMetadata(c(1,2,3),matrix(c(vstrMetadataOne, vstrMetadataTwo, vstrMetadataThree),ncol=length(vstrMetadataThree),byrow=TRUE),TRUE)

strLevelSelected1 = strsplit(lResult1[["names"]],"_")[[1]]
strLevelSelected1 = strLevelSelected1[length(strLevelSelected1)]
strLevelSelected2 = strsplit(lResult2[["names"]],"_")[[1]]
strLevelSelected2 = strLevelSelected2[length(strLevelSelected2)]
strLevelSelected3 = unlist(lapply(strsplit(lResult3[["names"]],"_"), function(x) x[[length(x)]]))
strLevelSelected4 = unlist(lapply(strsplit(lResult4[["names"]],"_"), function(x) x[[length(x)]]))

vxMetadata1 = matrix(rep(1, length(vstrMetadataOne)),nrow=1)
vxMetadata2 = matrix(rep(1, length(vstrMetadataOne)),nrow=1)

vxMetadata1[which(vstrMetadataOne==strLevelSelected1)] = 2
vstrNames1 = c(paste(c_strMetadata,3,"_",c_strLevel,"_",strLevelSelected1,sep=""))
vxMetadata2[which(vstrMetadataOne==strLevelSelected2)] = 2
vstrNames2 = c(paste(c_strMetadata,8,"_",c_strLevel,"_",strLevelSelected2,sep=""))

vxMetadata3 = matrix(c(vstrMetadataOne, vstrMetadataTwo, vstrMetadataThree),ncol=length(vstrMetadataThree),byrow=TRUE)
iLevelIndices = which(vxMetadata3==strLevelSelected3)
liDim = dim(vxMetadata3)
vxMetadata3[setdiff(c(1:(liDim[1]*liDim[2])),iLevelIndices)] = 1
vxMetadata3[iLevelIndices] = 2
vstrNames3 = c(paste(c_strMetadata,3,"_",c_strLevel,"_", strLevelSelected3[[1]],sep=""),paste(c_strMetadata,13,"_",c_strLevel,"_", strLevelSelected3[[2]],sep=""),paste(c_strMetadata,23,"_",c_strLevel,"_", strLevelSelected3[[3]],sep=""))

vxMetadata4 = matrix(c(vstrMetadataOne, vstrMetadataTwo, vstrMetadataThree),ncol=length(vstrMetadataThree),byrow=TRUE)
iLevelIndices = which(vxMetadata4==strLevelSelected4)
liDim = dim(vxMetadata4)
vxMetadata4[setdiff(c(1:(liDim[1]*liDim[2])),iLevelIndices)] = 1
vxMetadata4[iLevelIndices] = 2
vstrNames4 = c(paste(c_strMetadata,1,"_",c_strLevel,"_", strLevelSelected4[[1]],sep=""),paste(c_strMetadata,2,"_",c_strLevel,"_",strLevelSelected4[[2]],sep=""),paste(c_strMetadata,3,"_",c_strLevel,"_",strLevelSelected4[[3]],sep=""))

test_that("Test that factor data of different varying levels are correctly dummied and their names with levels are recoded.",{
  expect_equal(sort(lResult1[["names"]]), sort(vstrNames1))
  expect_equal(lResult1[["metadata"]], vxMetadata1)
  expect_equal(sort(lResult2[["names"]]), sort(vstrNames2))
  expect_equal(lResult2[["metadata"]], vxMetadata2)
  expect_equal(sort(lResult3[["names"]]), sort(vstrNames3))
  expect_equal(lResult3[["metadata"]], vxMetadata3)
  expect_equal(sort(lResult4[["names"]]), sort(vstrNames4))
  expect_equal(lResult4[["metadata"]], vxMetadata4)
})

context("Test funcPrepareMetadata with OUT dummying")
vstrMetadataOne = c(1,2,2,1,1,1,2,2,1,2,1,2)
vstrMetadataTwo = c(1,2,3,1,1,3,2,3,1,2,1,2)
vstrMetadataThree = c(1,2,3,4,1,3,2,4,1,2,1,2)

lResult1 = funcPrepareMetadata(c(3),vstrMetadataOne,FALSE)
lResult2 = funcPrepareMetadata(c(8),vstrMetadataOne,FALSE)
lResult3 = funcPrepareMetadata(c(3,13,23),matrix(c(vstrMetadataOne, vstrMetadataTwo, vstrMetadataThree),ncol=length(vstrMetadataThree),byrow=TRUE),FALSE)
lResult4 = funcPrepareMetadata(c(1,2,3),matrix(c(vstrMetadataOne, vstrMetadataTwo, vstrMetadataThree),ncol=length(vstrMetadataThree),byrow=TRUE),FALSE)

vstrNames1 = paste(c_strMetadata,c(3),sep="")
vstrNames2 = paste(c_strMetadata,c(8),sep="")
vstrNames3 = paste(c_strMetadata,c(3,13,23),sep="")
vstrNames4 = paste(c_strMetadata,c(1,2,3),sep="")

vxMetadata1 = vstrMetadataOne
vxMetadata2 = vxMetadata1
vxMetadata3 = matrix(c(vstrMetadataOne, vstrMetadataTwo, vstrMetadataThree),ncol=length(vstrMetadataThree),byrow=TRUE)
vxMetadata4 = matrix(c(vstrMetadataOne, vstrMetadataTwo, vstrMetadataThree),ncol=length(vstrMetadataThree),byrow=TRUE)

test_that("Test that factor data of different varying levels are not dummied and their names with levels are not recorded.",{
  expect_equal(lResult1[["names"]], vstrNames1)
  expect_equal(lResult1[["metadata"]], vxMetadata1)
  expect_equal(lResult2[["names"]], vstrNames2)
  expect_equal(lResult2[["metadata"]], vxMetadata2)
  expect_equal(lResult3[["names"]], vstrNames3)
  expect_equal(lResult3[["metadata"]], vxMetadata3)
  expect_equal(lResult4[["names"]], vstrNames4)
  expect_equal(lResult4[["metadata"]], vxMetadata4)
})


### funcQCSpikin
# PASS 10-22-2013
context("Test funcQCSpikin")
viMetadata1 = c(2,1,2,2,1,1,2,2,1,1)
viMetadata2 = c(2,1,1,1,2,2,1,1,1,1)
vdBug1 = c(0.1,0.3,0.2,0.6,0.8,0.9,0.3,0.5,0.6,0.4)
vdBug2 = c(0.1,0.0,0.2,0.0,0.8,0.0,0.3,0.0,0.6,0.0)
vdBug3 = c(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0)
test_that("Test that dummied data are correctly checked.",{
  expect_equal(funcQCSpikin(viMetadata1, vdBug1, 5)$PASS, TRUE)
  expect_equal(funcQCSpikin(viMetadata1, vdBug1, 3)$PASS, TRUE)
  expect_equal(funcQCSpikin(viMetadata1, vdBug1, 0)$PASS, TRUE)
  expect_equal(funcQCSpikin(viMetadata1, vdBug1, 6)$PASS, FALSE)
  expect_equal(funcQCSpikin(viMetadata1, vdBug1, 130)$PASS, FALSE)
  expect_equal(funcQCSpikin(viMetadata1, vdBug2, 2)$PASS, TRUE)
  expect_equal(funcQCSpikin(viMetadata1, vdBug2, 4)$PASS, FALSE)
  expect_equal(funcQCSpikin(viMetadata1, vdBug3, 1)$PASS, FALSE)
  expect_equal(funcQCSpikin(viMetadata2, vdBug1, 3)$PASS, TRUE)
  expect_equal(funcQCSpikin(viMetadata2, vdBug1, 2)$PASS, TRUE)
  expect_equal(funcQCSpikin(viMetadata2, vdBug1, 0)$PASS, TRUE)
  expect_equal(funcQCSpikin(viMetadata2, vdBug1, 4)$PASS, FALSE)
  expect_equal(funcQCSpikin(viMetadata2, vdBug1, 100)$PASS, FALSE)
  expect_equal(funcQCSpikin(viMetadata2, vdBug2, 2)$PASS, TRUE)
  expect_equal(funcQCSpikin(viMetadata2, vdBug2, 3)$PASS, FALSE)
  expect_equal(funcQCSpikin(viMetadata2, vdBug3, 1)$PASS, FALSE)
})

viMetadata3 = c(2,1,2,2,2,2,2,2,1,2)
viMetadata4 = rep(1,length(viMetadata1))
mtrxMatrixOf2 = matrix(c(viMetadata1,viMetadata2), byrow=TRUE, nrow=2)
mtrxMatrixOf3 = matrix(c(viMetadata3, viMetadata1,viMetadata2), byrow=TRUE, nrow=3)
mtrxMatrixOf4 = matrix(c(viMetadata3, viMetadata4, viMetadata1,viMetadata2), byrow=TRUE, nrow=4)
test_that("Test that dummied data are correctly check when they are matrices (row major).",{
  expect_equal(funcQCSpikin(mtrxMatrixOf2, vdBug1,3)$PASS, TRUE)
  expect_equal(funcQCSpikin(mtrxMatrixOf2, vdBug1,0)$PASS, TRUE)
  expect_equal(funcQCSpikin(mtrxMatrixOf2, vdBug1,4)$PASS, FALSE)
  expect_equal(funcQCSpikin(mtrxMatrixOf2, vdBug1,6)$PASS, FALSE)
  expect_equal(funcQCSpikin(mtrxMatrixOf3, vdBug1,3)$PASS, TRUE)
  expect_equal(funcQCSpikin(mtrxMatrixOf3, vdBug1,0)$PASS, TRUE)
  expect_equal(funcQCSpikin(mtrxMatrixOf3, vdBug1,4)$PASS, FALSE)
  expect_equal(funcQCSpikin(mtrxMatrixOf3, vdBug1,6)$PASS, FALSE)
  expect_equal(funcQCSpikin(mtrxMatrixOf4, vdBug1,3)$PASS, FALSE)
  expect_equal(funcQCSpikin(mtrxMatrixOf4, vdBug1,0)$PASS, TRUE)
  expect_equal(funcQCSpikin(mtrxMatrixOf4, vdBug1,4)$PASS, FALSE)
  expect_equal(funcQCSpikin(mtrxMatrixOf4, vdBug1,6)$PASS, FALSE)
  expect_equal(funcQCSpikin(mtrxMatrixOf2, vdBug2,1)$PASS, TRUE)
  expect_equal(funcQCSpikin(mtrxMatrixOf2, vdBug2,0)$PASS, TRUE)
  expect_equal(funcQCSpikin(mtrxMatrixOf2, vdBug2,4)$PASS, FALSE)
  expect_equal(funcQCSpikin(mtrxMatrixOf2, vdBug2,6)$PASS, FALSE)
  expect_equal(funcQCSpikin(mtrxMatrixOf3, vdBug2,2)$PASS, TRUE)
  expect_equal(funcQCSpikin(mtrxMatrixOf3, vdBug2,0)$PASS, TRUE)
  expect_equal(funcQCSpikin(mtrxMatrixOf3, vdBug2,4)$PASS, FALSE)
  expect_equal(funcQCSpikin(mtrxMatrixOf3, vdBug2,6)$PASS, FALSE)
  expect_equal(funcQCSpikin(mtrxMatrixOf4, vdBug2,3)$PASS, FALSE)
  expect_equal(funcQCSpikin(mtrxMatrixOf4, vdBug2,0)$PASS, TRUE)
  expect_equal(funcQCSpikin(mtrxMatrixOf4, vdBug2,4)$PASS, FALSE)
  expect_equal(funcQCSpikin(mtrxMatrixOf4, vdBug2,6)$PASS, FALSE)
})

viMetadata1 = c(1,2,1,1,2,2,1,1,2,2)
viMetadata2 = c(1,2,3,3,1,1,2,3,3,2)
vdBug1 = c(0.1,0.3,0.2,0.6,0.8,0.9,0.3,0.5,0.6,0.4)
vdBug2 = c(0.1,0.0,0.2,0.0,0.8,0.0,0.3,0.0,0.6,0.0)
vdBug3 = c(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0)

test_that("Test that dummied data are correctly checked without dummying.",{
  expect_equal(funcQCSpikin(viMetadata1, vdBug1, 5, FALSE)$PASS, TRUE)
  expect_equal(funcQCSpikin(viMetadata1, vdBug1, 3, FALSE)$PASS, TRUE)
  expect_equal(funcQCSpikin(viMetadata1, vdBug1, 0, FALSE)$PASS, TRUE)
  expect_equal(funcQCSpikin(viMetadata1, vdBug1, 12, FALSE)$PASS, FALSE)
  expect_equal(funcQCSpikin(viMetadata1, vdBug1, 130, FALSE)$PASS, FALSE)
  expect_equal(funcQCSpikin(viMetadata1, vdBug2, 2, FALSE)$PASS, TRUE)
  expect_equal(funcQCSpikin(viMetadata1, vdBug2, 6, FALSE)$PASS, FALSE)
  expect_equal(funcQCSpikin(viMetadata1, vdBug3, 1, FALSE)$PASS, FALSE)
  expect_equal(funcQCSpikin(viMetadata2, vdBug1, 9, FALSE)$PASS, TRUE)
  expect_equal(funcQCSpikin(viMetadata2, vdBug1, 2, FALSE)$PASS, TRUE)
  expect_equal(funcQCSpikin(viMetadata2, vdBug1, 0, FALSE)$PASS, TRUE)
  expect_equal(funcQCSpikin(viMetadata2, vdBug1, 12, FALSE)$PASS, FALSE)
  expect_equal(funcQCSpikin(viMetadata2, vdBug1, 100, FALSE)$PASS, FALSE)
  expect_equal(funcQCSpikin(viMetadata2, vdBug2, 3, FALSE)$PASS, TRUE)
  expect_equal(funcQCSpikin(viMetadata2, vdBug2, 9, FALSE)$PASS, FALSE)
  expect_equal(funcQCSpikin(viMetadata2, vdBug3, 1, FALSE)$PASS, FALSE)
})

viMetadata3 = c(1,1,2,2,3,3,4,4,5,5)
viMetadata4 = rep(1,length(viMetadata1))
mtrxMatrixOf2 = matrix(c(viMetadata1,viMetadata2), byrow=TRUE, nrow=2)
mtrxMatrixOf3 = matrix(c(viMetadata3, viMetadata1,viMetadata2), byrow=TRUE, nrow=3)
mtrxMatrixOf4 = matrix(c(viMetadata3, viMetadata4, viMetadata1,viMetadata2), byrow=TRUE, nrow=4)
test_that("Test that dummied data are correctly checked without dummying for matrices of data.",
{
  expect_equal(funcQCSpikin(mtrxMatrixOf2, vdBug1, 6, FALSE)$PASS, TRUE)
  expect_equal(funcQCSpikin(mtrxMatrixOf2, vdBug1, 12, FALSE)$PASS, FALSE)
  expect_equal(funcQCSpikin(mtrxMatrixOf2, vdBug1, 0, FALSE)$PASS, TRUE)
  expect_equal(funcQCSpikin(mtrxMatrixOf2, vdBug1, 100, FALSE)$PASS, FALSE)
  expect_equal(funcQCSpikin(mtrxMatrixOf3, vdBug1, 9, FALSE)$PASS, TRUE)
  expect_equal(funcQCSpikin(mtrxMatrixOf3, vdBug1, 8, FALSE)$PASS, TRUE)
  expect_equal(funcQCSpikin(mtrxMatrixOf3, vdBug1, 100, FALSE)$PASS, FALSE)
  expect_equal(funcQCSpikin(mtrxMatrixOf3, vdBug1, 0, FALSE)$PASS, TRUE)
  expect_equal(funcQCSpikin(mtrxMatrixOf4, vdBug1, 2, FALSE)$PASS, TRUE)
  expect_equal(funcQCSpikin(mtrxMatrixOf4, vdBug1, 8, FALSE)$PASS, TRUE)
  expect_equal(funcQCSpikin(mtrxMatrixOf4, vdBug1, 100, FALSE)$PASS, FALSE)
  expect_equal(funcQCSpikin(mtrxMatrixOf4, vdBug1, 0, FALSE)$PASS, TRUE)
  expect_equal(funcQCSpikin(mtrxMatrixOf2, vdBug2, 3, FALSE)$PASS, TRUE)
  expect_equal(funcQCSpikin(mtrxMatrixOf2, vdBug2, 12, FALSE)$PASS, FALSE)
  expect_equal(funcQCSpikin(mtrxMatrixOf2, vdBug2, 0, FALSE)$PASS, TRUE)
  expect_equal(funcQCSpikin(mtrxMatrixOf2, vdBug2, 100, FALSE)$PASS, FALSE)
  expect_equal(funcQCSpikin(mtrxMatrixOf3, vdBug2, 100, FALSE)$PASS, FALSE)
  expect_equal(funcQCSpikin(mtrxMatrixOf3, vdBug2, 3, FALSE)$PASS, TRUE)
  expect_equal(funcQCSpikin(mtrxMatrixOf3, vdBug2, 0, FALSE)$PASS, TRUE)
  expect_equal(funcQCSpikin(mtrxMatrixOf4, vdBug2, 100, FALSE)$PASS, FALSE)
  expect_equal(funcQCSpikin(mtrxMatrixOf4, vdBug2, 3, FALSE)$PASS, TRUE)
  expect_equal(funcQCSpikin(mtrxMatrixOf4, vdBug2, 0, FALSE)$PASS, TRUE)
})


### funcRoundMatrix
# PASS 10-22-2013
context("Test funcRoundMatrix for normal matrix and sparse matrix")
mtrxOne = matrix(c(1.1,2.2,3.3,4.4,5.5,0.6,7.7,8.8,9.9,10.10,11.11,0.12,13.13,14.14,15.15,16.16,0.17,18.18,19.19,20.20,21.21),nrow=7,ncol=3)
mtrxTwo = matrix(c(1.1,0,3.3,4.4,0,0.6,7.7,8.8,9.9,0,11.11,0.12,13.13,14.14,0,16.16,0.17,18.18,19.19,20.20,21.21),nrow=7,ncol=3)
# Not treated like zero inflated
mtrxOneAnswer = matrix(c(1,2,3,4,6,1,8,9,10,10,11,0,13,14,15,16,0,18,19,20,21),nrow=7,ncol=3)
mtrxTwoAnswer = matrix(c(1,0,3,4,0,1,8,9,10,0,11,0,13,14,0,16,0,18,19,20,21),nrow=7,ncol=3)
# Treated like zero inflated
mtrxOneAnswer2 = matrix(c(1,2,3,4,6,1,8,9,10,10,11,1,13,14,15,16,1,18,19,20,21),nrow=7,ncol=3)
mtrxTwoAnswer2 = matrix(c(1,0,3,4,0,1,8,9,10,0,11,1,13,14,0,16,1,18,19,20,21),nrow=7,ncol=3)

test_that("Test that rounding happens as normal expect no 0s are added, 1 is the lowest number rounding can occur to.",
{
  expect_equal( mtrxOneAnswer, funcRoundMatrix( mtrxOne, FALSE ) )
  expect_equal( mtrxOneAnswer2, funcRoundMatrix( mtrxOne, TRUE ) )
  expect_equal( mtrxTwoAnswer, funcRoundMatrix( mtrxTwo, FALSE ) )
  expect_equal( mtrxTwoAnswer2, funcRoundMatrix( mtrxTwo, TRUE ) )
})


### funcSample
# PASS 10-22-2013
context( "Test to make sure behavior of vectorize sample is consistent between 1 element and multiple element vectors." )
vOne = 1:10
iLengthOne = length(vOne)

test_that( "Test that vector of multiple samples sample correctly.",
{
  expect_equal( vOne, sort(funcSample(vOne, iLengthOne, replace = FALSE ) ) )
  expect_equal( iLengthOne/2, length( sort( funcSample( vOne, iLengthOne/2, replace = TRUE ) ) ) )
  expect_equal( iLengthOne*2, length( sort( funcSample( vOne, iLengthOne*2, replace = TRUE ) ) ) )
  expect_true( unique( funcSample( vOne, iLengthOne*2, replace = TRUE ) %in% vOne ) )
})

test_that( "Test that 1 element is sampled correctly.",
{
  expect_equal( 1, funcSample( c( 1 ), 1, replace = FALSE ) )
  expect_equal( rep( 1, 3 ), funcSample( c( 1 ), 3, replace = TRUE ) )
  expect_equal( 11, funcSample( c( 11 ), 1, replace = FALSE ) )
  expect_equal( rep( 11, 3 ), funcSample( c( 11 ), 3, replace = TRUE ) )
})

### funcShuffleMatrix
# PASS 10-22-2013
context( "Test funcShuffleMatrix for even shuffling and read depth improvement" )
mtrxOne = matrix(c(1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0),nrow=7,ncol=3)
test_that("Test that the matrices shuffle evenly and leave no samples of sum 0.",
{
  expect_true( unique( colSums( funcShuffleMatrix( mtrxOne, 1 ) ) %in% c( 5, 1 ) ) )
  expect_true( unique( colSums( funcShuffleMatrix( mtrxOne, 2 ) ) %in% c( 2, 3 ) ) )
  expect_true( unique( colSums( funcShuffleMatrix( mtrxOne, 7 ) ) %in% c( 7, 0 ) ) )
})


### funcSpikeNewBug
# PASS 10-22-2013
# funcSpikeNewBug( vdContinuousMetadata1, vdNotSparseBug1, iMultiplier1 )
context("Test funcSpikeNewBug")

test_that("funcSpikeNewBug: Test that relationships are made with simple inputs.",{
  expect_equal( 1:10, funcSpikeNewBug( 1:10, 1:10, 1 ))
  expect_equal( 1:10, funcSpikeNewBug( 1:10, 1:10, 5 ))
  expect_equal( rep(1,8), funcSpikeNewBug( rep(1,8), rep(1,8), 1 ))
  expect_equal( rep(1,8), funcSpikeNewBug( rep(1,8), rep(1,8), 5 ))
  expect_equal( rep(5.5,10), funcSpikeNewBug( 1:10, 10:1, 1 ))
  expect_equal( c(9.1, 8.3, 7.5, 6.7, 5.9, 5.1, 4.3, 3.5, 2.7, 1.9), funcSpikeNewBug( 10:1, 1:10, 9 ))
  expect_equal( 1:10, funcSpikeNewBug( matrix(c(1:10,1:10), byrow=TRUE,nrow=2,ncol=10), 1:10, 1 ))
  expect_equal( 1:10, funcSpikeNewBug( matrix(c(1:10,1:10), byrow=TRUE,nrow=2,ncol=10), 1:10, 5 ))
  expect_equal( round(c(2.285714, 2.428571, 2.571429, 2.714286),2), round(funcSpikeNewBug( matrix(c(1:4,4:1), byrow=TRUE,nrow=2,ncol=4), 1:4, 3 ),2))
  expect_equal( which(c(1,0,3,4,5,6,0,8,9,0)==0), which(funcSpikeNewBug( 1:10, c(1,0,3,4,5,6,0,8,9,0),1)==0))
  expect_equal( which(c(1,0,1,1,0,1,0,1,0,0)==0), which(funcSpikeNewBug( matrix(c(1:10,1:10), byrow=TRUE,nrow=2,ncol=10), c(10,0,20,10,0,10,0,10,0,0), 5 )==0))
})


### funcTruncatedRLNorm
# PASS 10-22-2013
iCount <- 1000000
dMu <- 2
dSD <- 3
dExp = funcGetExp( 2, 3 )
iThreshold <- NA
lResult <- funcTruncatedRLNorm( iCount, log( dMu ), log( dSD ), iThreshold )
test_that( "Test that funcTruncatedRLNorm creates an expected distribution without truncation.", {
  expect_true( round( dExp ) - 1 < round( mean( lResult ) ) && round( dExp ) + 1 > round( mean( lResult ) ) )
  expect_equal( iCount, length( lResult ) )
})

iCount <- 1000000
dMu <- 5
dSD <- 7
dExp = funcGetExp( 5, 7 )
iThreshold <- NA
lResult <- funcTruncatedRLNorm( iCount, log( dMu ), log( dSD ), iThreshold )
test_that( "Test that funcTruncatedRLNorm creates an expected distribution without truncation.", {
  expect_true( round( dExp ) - 1 < round( mean( lResult ) ) && round( dExp ) + 1 > round( mean( lResult ) ) )
  expect_equal( iCount, length( lResult ) )
})

iCount <- 100
dMu <- 2
dSD <- 3
dExp = funcGetExp( 2, 3 )
iThreshold <- 3
lResult <- funcTruncatedRLNorm( iCount, log( dMu ), log( dSD ), iThreshold )
test_that( "Test that funcTruncatedRLNorm creates an expected distribution without truncation.", {
  expect_equal( 0, length( lResult[ lResult > iThreshold ] ) )
})

iCount <- 1000
dMu <- 2
dSD <- 7
dExp = funcGetExp( 2, 7 )
iThreshold <- 3
lResult <- funcTruncatedRLNorm( iCount, log( dMu ), log( dSD ), iThreshold )
test_that( "Test that funcTruncatedRLNorm creates an expected distribution without truncation." ,{
  expect_equal( 0, length( lResult[ lResult > iThreshold ] ) )
})


### funcUpdateDistributionToExpectation
# Pass 10-22-2013
vdDist1 = rlnorm( 10, log( 2 ), log( 3 ) )
dMean1 = round( mean( vdDist1 ) )
dTest1 = round( mean( funcUpdateDistributionToExpectation( vdDist1, floor( dMean1 + 10 ) ) ) )
# Tests that minimum expectation is 1
dTest2 = round( mean( funcUpdateDistributionToExpectation( vdDist1, floor( dMean1 - 10 ) ) ) )
vdDist2 = rlnorm( 100, log( 10 ), log( 3 ) )
dMean2 = round( mean( vdDist2 ) )
dTest3 = round( mean( funcUpdateDistributionToExpectation( vdDist2, floor( dMean2 + 10 ) ) ) )
dTest4 = round( mean( funcUpdateDistributionToExpectation( vdDist2, floor( dMean2 - 10 ) ) ) )
test_that( "Test that the funcUpdateDistributionToExpectation, updates distributions to the correct mean.", {
  expect_true( ( dTest1 >= floor( dMean1 + 10 ) - 1 ) && ( dTest1 <= floor( dMean1 + 10 ) + 1 ) )
  expect_true( dTest2 == 1 )
  expect_true( ( dTest3 >= floor( dMean2 + 10 ) - 1 ) && ( dTest3 <= floor( dMean2 + 10 ) + 1 ) )
  expect_true( ( dTest4 >= floor( dMean2 - 10 ) - 1 ) && ( dTest4 <= floor( dMean2 - 10 ) + 1 ) )
})


### funcUpdateDistributionToSum
# PASS 10-22-2013
vdDist1 = rlnorm( 10, log( 2 ), log( 3 ) )
dSum1 = round( sum( vdDist1 ) )
dTest1 = round( sum( funcUpdateDistributionToSum( vdDist1, floor( dSum1 + 10 ) ) ) )
dTest2 = round( sum( funcUpdateDistributionToSum( vdDist1, floor( dSum1 - 10 ) ) ) )
vdDist2 = rlnorm( 100, log( 2 ), log( 3 ) )
dSum2 = round( sum( vdDist2 ) )
dTest3 = round( sum( funcUpdateDistributionToSum( vdDist2, floor( dSum2 + 10 ) ) ) )
dTest4 = round( sum( funcUpdateDistributionToSum( vdDist2, floor( dSum2 - 10 ) ) ) )
test_that( "Test that the funcUpdateDistributionToSum, updates distributions to the correct sum.", {
  expect_true( ( dTest1 >= floor( dSum1 + 10 ) - 1 ) && ( dTest1 <= floor( dSum1 + 10 ) + 1 ) )
  expect_true( ( dTest2 >= floor( dSum1 - 10 ) - 1 ) && ( dTest2 <= floor( dSum1 - 10 ) + 1 ) )
  expect_true( ( dTest3 >= floor( dSum2 + 10 ) - 1 ) && ( dTest3 <= floor( dSum2 + 10 ) + 1 ) )
  expect_true( ( dTest4 >= floor( dSum2 - 10 ) - 1 ) && ( dTest4 <= floor( dSum2 - 10 ) + 1 ) )
})


### func_zero_inflate
# PASS 10-22-2013
context("Test func_zero_inflate checking that the made distributions are equal to the parameters.")
# No SD
dMeanOne = 5
dPercentZeroOne = .50
dSDOne = 1
iNumberSamples = 100000
vdTestOne = func_zero_inflate(dMeanOne,dPercentZeroOne,iNumberSamples,dSDOne)
viNotZeroOne = which(vdTestOne!=0)
dRealMeanOne = mean(log(vdTestOne[viNotZeroOne]))
dRealSDOne = sd(log(vdTestOne[viNotZeroOne]))
iRealLengthOne = length(vdTestOne)
dRealPercentZeroOne = length(which(vdTestOne==0))/iRealLengthOne
# SD
dMeanTwo = 10
dPercentZeroTwo = .3
dSDTwo = 2
iNumberSamples = 100000
vdTestTwo = func_zero_inflate(dMeanTwo,dPercentZeroTwo,iNumberSamples,dSDTwo)
viNotZeroTwo = which(vdTestTwo!=0)
dRealMeanTwo = mean(log(vdTestTwo[viNotZeroTwo]))
dRealSDTwo = sd(log(vdTestTwo[viNotZeroTwo]))
iRealLengthTwo = length(vdTestTwo)
dRealPercentZeroTwo = length(which(vdTestTwo==0))/iRealLengthTwo
# No negatives
dMeanThree = 1
dPercentZeroThree = .0
dSDThree = 5
iNumberSamples = 100000
vdTestThree = func_zero_inflate(dMeanThree,dPercentZeroThree,iNumberSamples,dSDThree)

test_that("Test that func_zero_inflate creates a distribution close to the given parameters.",{
  expect_equal(round(dRealMeanOne),round(dMeanOne))
  expect_equal(round(dRealSDOne),round(dSDOne))
  expect_equal(iRealLengthOne,iNumberSamples)
  expect_equal(round(dPercentZeroOne,1),round(dRealPercentZeroOne,1))

  expect_equal(round(dRealMeanTwo),round(dMeanTwo))
  expect_equal(round(dRealSDTwo),round(dSDTwo))
  expect_equal(iRealLengthTwo,iNumberSamples)
  expect_equal(round(dPercentZeroTwo,1),round(dRealPercentZeroTwo,1))

  expect_equal(length(which(vdTestThree<0)),0)
})