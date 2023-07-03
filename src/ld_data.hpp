/*  _______________________________________________________________________

    DAKOTA: Design Analysis Kit for Optimization and Terascale Applications
    Copyright 2014-2023
    National Technology & Engineering Solutions of Sandia, LLC (NTESS).
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Dakota directory.
    _______________________________________________________________________ */

//- Class:	     
//- Description: Data for low-discrepancy generators
//- Owner:       Pieterjan Robbe
//- Checked by:
//- Version:

#ifndef LD_DATA_H
#define LD_DATA_H

#include "dakota_data_types.hpp"

namespace Dakota {

/// Generating vector from Cools, Kuo and Nuyens (2006)
UInt32 cools_kuo_nuyens_d250_m20[] = {
  1,
  182667,
  469891,
  498753,
  110745,
  446247,
  250185,
  118627,
  245333,
  283199,
  408519,
  391023,
  246327,
  126539,
  399185,
  461527,
  300343,
  69681,
  516695,
  436179,
  106383,
  238523,
  413283,
  70841,
  47719,
  300129,
  113029,
  123925,
  410745,
  211325,
  17489,
  511893,
  40767,
  186077,
  519471,
  255369,
  101819,
  243573,
  66189,
  152143,
  503455,
  113217,
  132603,
  463967,
  297717,
  157383,
  224015,
  502917,
  36237,
  94049,
  170665,
  79397,
  123963,
  223451,
  323871,
  303633,
  98567,
  318855,
  494245,
  477137,
  177975,
  64483,
  26695,
  88779,
  94497,
  239429,
  381007,
  110205,
  339157,
  73397,
  407559,
  181791,
  442675,
  301397,
  32569,
  147737,
  189949,
  138655,
  350241,
  63371,
  511925,
  515861,
  434045,
  383435,
  249187,
  492723,
  479195,
  84589,
  99703,
  239831,
  269423,
  182241,
  61063,
  130789,
  143095,
  471209,
  139019,
  172565,
  487045,
  304803,
  45669,
  380427,
  19547,
  425593,
  337729,
  237863,
  428453,
  291699,
  238587,
  110653,
  196113,
  465711,
  141583,
  224183,
  266671,
  169063,
  317617,
  68143,
  291637,
  263355,
  427191,
  200211,
  365773,
  254701,
  368663,
  248047,
  209221,
  279201,
  323179,
  80217,
  122791,
  316633,
  118515,
  14253,
  129509,
  410941,
  402601,
  511437,
  10469,
  366469,
  463959,
  442841,
  54641,
  44167,
  19703,
  209585,
  69037,
  33317,
  433373,
  55879,
  245295,
  10905,
  468881,
  128617,
  417919,
  45067,
  442243,
  359529,
  51109,
  290275,
  168691,
  212061,
  217775,
  405485,
  313395,
  256763,
  152537,
  326437,
  332981,
  406755,
  423147,
  412621,
  362019,
  279679,
  169189,
  107405,
  251851,
  5413,
  316095,
  247945,
  422489,
  2555,
  282267,
  121027,
  369319,
  204587,
  445191,
  337315,
  322505,
  388411,
  102961,
  506099,
  399801,
  254381,
  452545,
  309001,
  147013,
  507865,
  32283,
  320511,
  264647,
  417965,
  227069,
  341461,
  466581,
  386241,
  494585,
  201479,
  151243,
  481337,
  68195,
  75401,
  58359,
  448107,
  459499,
  9873,
  365117,
  350845,
  181873,
  7917,
  436695,
  43899,
  348367,
  423927,
  437399,
  385089,
  21693,
  268793,
  49257,
  250211,
  125071,
  341631,
  310163,
  94631,
  108795,
  21175,
  142847,
  383599,
  71105,
  65989,
  446433,
  177457,
  107311,
  295679,
  442763,
  40729,
  322721,
  420175,
  430359,
  480757
};

/// Recommended generating vector from Kuo's website
/// \url{https://web.maths.unsw.edu.au/~fkuo/lattice/index.html}
UInt32 kuo_d3600_m20[] = {
  1,
  182667,
  469891,
  498753,
  110745,
  446247,
  250185,
  118627,
  245333,
  283199,
  408519,
  391023,
  246327,
  126539,
  399185,
  461527,
  300343,
  69681,
  516695,
  436179,
  106383,
  238523,
  413283,
  70841,
  47719,
  300129,
  113029,
  123925,
  410745,
  211325,
  17489,
  511893,
  40767,
  186077,
  519471,
  255369,
  101819,
  243573,
  66189,
  152143,
  503455,
  113217,
  132603,
  463967,
  297717,
  157383,
  224015,
  502917,
  36237,
  94049,
  170665,
  79397,
  123963,
  223451,
  323871,
  303633,
  98567,
  318855,
  494245,
  477137,
  177975,
  64483,
  26695,
  88779,
  94497,
  239429,
  381007,
  110205,
  339157,
  73397,
  407559,
  181791,
  442675,
  301397,
  32569,
  147737,
  189949,
  138655,
  350241,
  63371,
  511925,
  515861,
  434045,
  383435,
  249187,
  492723,
  479195,
  84589,
  99703,
  239831,
  269423,
  182241,
  61063,
  130789,
  143095,
  471209,
  139019,
  172565,
  487045,
  304803,
  45669,
  380427,
  19547,
  425593,
  337729,
  237863,
  428453,
  291699,
  238587,
  110653,
  196113,
  465711,
  141583,
  224183,
  266671,
  169063,
  317617,
  68143,
  291637,
  263355,
  427191,
  200211,
  365773,
  254701,
  368663,
  248047,
  209221,
  279201,
  323179,
  80217,
  122791,
  316633,
  118515,
  14253,
  129509,
  410941,
  402601,
  511437,
  10469,
  366469,
  463959,
  442841,
  54641,
  44167,
  19703,
  209585,
  69037,
  33317,
  433373,
  55879,
  245295,
  10905,
  468881,
  128617,
  417919,
  45067,
  442243,
  359529,
  51109,
  290275,
  168691,
  212061,
  217775,
  405485,
  313395,
  256763,
  152537,
  326437,
  332981,
  406755,
  423147,
  412621,
  362019,
  279679,
  169189,
  107405,
  251851,
  5413,
  316095,
  247945,
  422489,
  2555,
  282267,
  121027,
  369319,
  204587,
  445191,
  337315,
  322505,
  388411,
  102961,
  506099,
  399801,
  254381,
  452545,
  309001,
  147013,
  507865,
  32283,
  320511,
  264647,
  417965,
  227069,
  341461,
  466581,
  386241,
  494585,
  201479,
  151243,
  481337,
  68195,
  75401,
  58359,
  448107,
  459499,
  9873,
  365117,
  350845,
  181873,
  7917,
  436695,
  43899,
  348367,
  423927,
  437399,
  385089,
  21693,
  268793,
  49257,
  250211,
  125071,
  341631,
  310163,
  94631,
  108795,
  21175,
  142847,
  383599,
  71105,
  65989,
  446433,
  177457,
  107311,
  295679,
  442763,
  40729,
  322721,
  420175,
  430359,
  480757,
  221417,
  460963,
  263529,
  263625,
  16919,
  289121,
  243609,
  414769,
  147195,
  326731,
  277925,
  185497,
  29735,
  11185,
  212811,
  371743,
  33547,
  250809,
  390965,
  486257,
  348223,
  263649,
  241785,
  425685,
  240127,
  447899,
  182195,
  306733,
  465107,
  434787,
  239361,
  437289,
  183495,
  75179,
  431157,
  376487,
  141253,
  293903,
  399797,
  325045,
  277541,
  215657,
  263115,
  193361,
  478441,
  299477,
  319217,
  176929,
  279451,
  32237,
  346181,
  374911,
  443279,
  319283,
  306211,
  70449,
  88145,
  272219,
  178705,
  159715,
  177793,
  164999,
  16833,
  49379,
  83863,
  223583,
  22615,
  225323,
  111767,
  361309,
  182521,
  154685,
  426543,
  482481,
  122871,
  250773,
  324031,
  227329,
  465853,
  165123,
  163615,
  245155,
  192403,
  364635,
  276575,
  195289,
  439767,
  386711,
  472741,
  6385,
  435121,
  299865,
  151695,
  344227,
  505437,
  130021,
  254671,
  165583,
  129277,
  513679,
  155197,
  147197,
  161911,
  408721,
  134261,
  56599,
  410699,
  110337,
  18283,
  393383,
  222965,
  244381,
  89807,
  198515,
  219643,
  467277,
  69251,
  227483,
  90213,
  405587,
  225651,
  349029,
  365859,
  448277,
  284285,
  85247,
  8363,
  235701,
  168505,
  11797,
  326593,
  482271,
  181411,
  424383,
  486111,
  136879,
  467187,
  523659,
  357769,
  102621,
  120121,
  253107,
  11107,
  410147,
  315347,
  1551,
  60525,
  265467,
  404483,
  174571,
  303643,
  12741,
  486785,
  382527,
  264167,
  438949,
  73635,
  149281,
  482205,
  150185,
  256827,
  108249,
  38907,
  511485,
  437613,
  401415,
  196935,
  228179,
  207375,
  366367,
  55589,
  84509,
  108335,
  371355,
  384565,
  35521,
  175225,
  111363,
  61861,
  183569,
  406605,
  71801,
  211703,
  132611,
  241903,
  215337,
  235709,
  415963,
  274405,
  4703,
  240525,
  416779,
  243755,
  157423,
  286867,
  227441,
  254477,
  462447,
  276941,
  111189,
  144789,
  368511,
  321515,
  355653,
  330037,
  193619,
  172623,
  357807,
  255105,
  265789,
  348989,
  399829,
  88709,
  147959,
  397901,
  25605,
  426447,
  223367,
  281913,
  109421,
  354737,
  97901,
  221801,
  483127,
  331811,
  284329,
  200029,
  472285,
  161795,
  142317,
  521093,
  375505,
  103755,
  166279,
  57527,
  161091,
  215319,
  485895,
  25887,
  89591,
  240303,
  253535,
  491933,
  137349,
  26775,
  1729,
  409481,
  468105,
  114849,
  396207,
  131517,
  179647,
  82499,
  270263,
  507511,
  126811,
  392021,
  290691,
  144999,
  409189,
  25557,
  217529,
  292351,
  418783,
  192003,
  162237,
  300457,
  283151,
  244869,
  372931,
  245703,
  124071,
  1135,
  286385,
  256437,
  277087,
  351965,
  88117,
  379293,
  177055,
  5879,
  372417,
  277241,
  204087,
  28287,
  510057,
  509695,
  148637,
  173693,
  18863,
  410971,
  490307,
  158579,
  180629,
  292825,
  56021,
  161863,
  391107,
  441385,
  41333,
  88349,
  248261,
  414605,
  250945,
  258519,
  519295,
  428331,
  82417,
  297583,
  236509,
  193353,
  521611,
  252473,
  8615,
  315077,
  356463,
  168977,
  225931,
  471707,
  523789,
  283671,
  107083,
  338195,
  64035,
  422417,
  114519,
  179967,
  228957,
  249611,
  300357,
  414135,
  316269,
  482363,
  349081,
  247739,
  257763,
  514731,
  459623,
  273057,
  190725,
  16075,
  225093,
  263631,
  304121,
  456955,
  185759,
  447357,
  463223,
  507381,
  91535,
  293531,
  336195,
  338591,
  233527,
  475567,
  202419,
  497965,
  427305,
  292097,
  87701,
  46525,
  425127,
  162119,
  91059,
  509853,
  439515,
  193197,
  489889,
  247123,
  72951,
  19091,
  32249,
  175185,
  401645,
  481103,
  156391,
  461987,
  274743,
  258789,
  363861,
  329821,
  441411,
  494409,
  523681,
  275133,
  476635,
  356313,
  301287,
  370055,
  223247,
  474385,
  235885,
  464331,
  100449,
  154789,
  354597,
  46075,
  509469,
  293887,
  118997,
  94605,
  130853,
  62381,
  446211,
  110775,
  441749,
  11647,
  441327,
  518019,
  266745,
  303321,
  428073,
  439511,
  327309,
  370317,
  99099,
  92823,
  47303,
  508365,
  122925,
  371521,
  413947,
  381997,
  128693,
  338147,
  68609,
  298361,
  142827,
  78791,
  110299,
  74983,
  399743,
  248095,
  446313,
  211653,
  42373,
  380463,
  305563,
  287069,
  174685,
  350981,
  87739,
  118153,
  130365,
  288001,
  85703,
  370809,
  261363,
  243289,
  127027,
  450185,
  230997,
  96425,
  165333,
  329901,
  166631,
  451049,
  124935,
  425907,
  187583,
  426405,
  407503,
  481683,
  436815,
  36701,
  230675,
  142295,
  384337,
  97869,
  174277,
  105907,
  411309,
  442025,
  420455,
  275357,
  36543,
  399931,
  369307,
  300069,
  198191,
  275635,
  74869,
  395981,
  154569,
  425805,
  209447,
  511999,
  311137,
  401091,
  114095,
  123237,
  409775,
  313845,
  411497,
  307259,
  252635,
  486853,
  415967,
  393585,
  252305,
  248807,
  429869,
  19877,
  238767,
  465765,
  251435,
  31931,
  9375,
  394037,
  482327,
  226155,
  276901,
  145515,
  515957,
  58865,
  432665,
  317515,
  202675,
  194919,
  363725,
  167451,
  496379,
  88081,
  429401,
  193909,
  464257,
  21783,
  255519,
  173419,
  4525,
  396301,
  323707,
  301965,
  58913,
  181921,
  392791,
  517181,
  199121,
  150613,
  498367,
  353209,
  362387,
  426477,
  30755,
  481743,
  79335,
  242545,
  193489,
  154217,
  152365,
  248003,
  160629,
  450905,
  119515,
  312019,
  221569,
  470791,
  42569,
  305501,
  107903,
  322337,
  40815,
  2239,
  45407,
  377109,
  222897,
  408051,
  334611,
  102593,
  474007,
  449755,
  395473,
  183207,
  445655,
  359145,
  480015,
  506445,
  415163,
  459765,
  240993,
  348643,
  311311,
  301277,
  506281,
  478109,
  106897,
  475519,
  295603,
  329847,
  177575,
  127243,
  272713,
  304281,
  50573,
  349417,
  131487,
  502907,
  412417,
  397071,
  147315,
  235781,
  255885,
  421135,
  161223,
  77175,
  331995,
  78655,
  211715,
  301259,
  219035,
  127701,
  83493,
  61371,
  300739,
  312513,
  466155,
  314441,
  345733,
  156889,
  487233,
  146325,
  73709,
  23791,
  469019,
  429905,
  210955,
  512205,
  237585,
  38469,
  380697,
  299875,
  275393,
  388021,
  354889,
  266469,
  293501,
  186931,
  124695,
  432785,
  242425,
  107199,
  57331,
  200461,
  371767,
  159645,
  217483,
  299463,
  443487,
  387823,
  288247,
  231339,
  69451,
  167207,
  474283,
  288285,
  35023,
  310071,
  52017,
  377933,
  184317,
  236333,
  100811,
  6867,
  225535,
  495413,
  237521,
  165557,
  466251,
  494583,
  41971,
  325869,
  497911,
  107613,
  261385,
  254175,
  420485,
  345119,
  271503,
  474889,
  79655,
  131139,
  494331,
  204103,
  56679,
  30945,
  103137,
  387471,
  375569,
  347443,
  84925,
  109079,
  121365,
  28371,
  23283,
  402765,
  45085,
  203239,
  332045,
  280521,
  511521,
  420779,
  457787,
  308723,
  213891,
  340515,
  6077,
  469349,
  317283,
  193051,
  84615,
  26411,
  111133,
  384309,
  89225,
  15061,
  451299,
  253109,
  283983,
  272139,
  224493,
  457497,
  198591,
  60857,
  170653,
  516063,
  63655,
  176635,
  496653,
  142945,
  181199,
  271031,
  262767,
  431349,
  82049,
  441275,
  496297,
  469933,
  413335,
  194583,
  105247,
  500683,
  139889,
  100381,
  440753,
  193069,
  5979,
  193247,
  301363,
  33843,
  221103,
  132691,
  130351,
  367705,
  132071,
  410455,
  5763,
  280417,
  359827,
  30087,
  290263,
  357573,
  86205,
  65473,
  144355,
  146417,
  84499,
  5649,
  518281,
  65049,
  42149,
  145689,
  314395,
  235127,
  205711,
  182221,
  346969,
  227929,
  428517,
  31343,
  431861,
  300643,
  498615,
  100851,
  488335,
  310185,
  322423,
  497201,
  484405,
  444527,
  51333,
  413769,
  65027,
  412357,
  287663,
  99387,
  508665,
  134905,
  337535,
  351489,
  290759,
  513485,
  391501,
  139009,
  514003,
  360051,
  194857,
  214625,
  190319,
  193121,
  279099,
  188459,
  465397,
  444741,
  1599,
  376567,
  502721,
  300331,
  309899,
  105713,
  221565,
  209529,
  226541,
  168083,
  172953,
  480521,
  41525,
  503357,
  513735,
  20585,
  509665,
  415835,
  268953,
  318039,
  318647,
  61943,
  406283,
  159249,
  243265,
  86427,
  356525,
  161945,
  88741,
  155739,
  30999,
  469509,
  416535,
  516925,
  20595,
  492481,
  118079,
  67531,
  93975,
  254945,
  460265,
  32577,
  286941,
  461315,
  229511,
  205789,
  123411,
  228447,
  405775,
  316867,
  12567,
  480181,
  507221,
  463329,
  360803,
  271225,
  433273,
  489915,
  164057,
  389401,
  61883,
  365689,
  120717,
  7725,
  97135,
  145163,
  311659,
  517775,
  47851,
  109593,
  472773,
  400637,
  19087,
  293147,
  484279,
  38049,
  331861,
  226935,
  451367,
  60083,
  507751,
  256185,
  159613,
  478559,
  466397,
  243859,
  474635,
  69159,
  75285,
  471867,
  436483,
  194003,
  332149,
  521247,
  221209,
  459463,
  195993,
  201455,
  15545,
  137185,
  204247,
  474959,
  324651,
  121031,
  328369,
  486033,
  104301,
  256295,
  306443,
  275897,
  321283,
  350205,
  330847,
  367481,
  475923,
  207701,
  494919,
  336813,
  267257,
  227219,
  472323,
  78973,
  513131,
  225923,
  434771,
  483219,
  114577,
  357429,
  165589,
  408211,
  494065,
  77643,
  355371,
  286569,
  35217,
  179467,
  36367,
  443261,
  142009,
  354987,
  186625,
  817,
  243731,
  383875,
  40125,
  13155,
  367409,
  468829,
  447387,
  13997,
  356469,
  285679,
  92327,
  419925,
  442095,
  392039,
  176961,
  146825,
  175823,
  88141,
  199553,
  14927,
  123441,
  90739,
  241217,
  53023,
  264725,
  283201,
  143723,
  408395,
  481677,
  348225,
  320161,
  37003,
  522423,
  70665,
  340285,
  146281,
  445675,
  385085,
  250727,
  502921,
  399081,
  199841,
  400747,
  180879,
  442621,
  338349,
  376305,
  286299,
  298983,
  148231,
  188293,
  312811,
  374571,
  15025,
  80623,
  460867,
  71815,
  191007,
  370419,
  285835,
  103895,
  34305,
  366809,
  137355,
  312313,
  393671,
  123371,
  211783,
  492741,
  371635,
  68305,
  448925,
  491667,
  390687,
  469303,
  455913,
  131369,
  23171,
  387865,
  229647,
  357387,
  44221,
  85497,
  497103,
  331905,
  63131,
  522595,
  121011,
  512597,
  367289,
  91267,
  72693,
  40989,
  83659,
  344031,
  323869,
  453437,
  226351,
  29011,
  274231,
  396789,
  439683,
  397207,
  341859,
  374821,
  213971,
  20097,
  174827,
  437687,
  205455,
  109395,
  309329,
  32869,
  284877,
  333319,
  167233,
  157959,
  389749,
  484567,
  265265,
  99497,
  288789,
  319989,
  133933,
  487021,
  278713,
  323331,
  450597,
  404371,
  135937,
  311647,
  363539,
  486493,
  440019,
  414487,
  232917,
  851,
  429267,
  93357,
  383357,
  112977,
  390981,
  139493,
  117513,
  482133,
  447549,
  14093,
  218813,
  320407,
  212417,
  34553,
  38101,
  185693,
  185261,
  474591,
  502403,
  204315,
  476157,
  118695,
  262023,
  217349,
  433597,
  276007,
  205027,
  497981,
  159879,
  450003,
  78073,
  449741,
  394723,
  276655,
  302909,
  73193,
  485803,
  128453,
  190853,
  448625,
  423923,
  450751,
  84341,
  88515,
  278945,
  15023,
  460179,
  432595,
  31905,
  81855,
  83211,
  460571,
  125237,
  2637,
  206799,
  382775,
  96453,
  467767,
  96071,
  332519,
  175269,
  187643,
  337227,
  378167,
  6565,
  317723,
  435053,
  39425,
  377629,
  240171,
  210049,
  141119,
  341063,
  318885,
  239455,
  110491,
  65875,
  335251,
  336349,
  77589,
  157821,
  62429,
  522025,
  168735,
  130331,
  489131,
  446475,
  116305,
  262001,
  99991,
  482961,
  421499,
  245947,
  359321,
  38481,
  321173,
  161151,
  406195,
  7507,
  72795,
  256199,
  306699,
  314195,
  438761,
  181227,
  349763,
  431907,
  26319,
  441141,
  44541,
  486865,
  326673,
  379261,
  361973,
  521935,
  19205,
  90667,
  231881,
  287905,
  269053,
  128239,
  251485,
  423591,
  503993,
  225627,
  155339,
  518847,
  374369,
  329871,
  219073,
  360565,
  54941,
  244749,
  414575,
  281531,
  229491,
  321177,
  287485,
  301621,
  206741,
  340997,
  231751,
  389515,
  347819,
  176747,
  161639,
  429653,
  445151,
  94467,
  224415,
  291995,
  497733,
  414631,
  418071,
  233485,
  470293,
  96223,
  76519,
  181327,
  105741,
  25159,
  491531,
  160453,
  314353,
  99767,
  337271,
  113963,
  467495,
  30361,
  499149,
  13275,
  67849,
  351547,
  452929,
  237379,
  392947,
  337257,
  190073,
  294111,
  201073,
  216495,
  25003,
  513811,
  316363,
  211719,
  258485,
  274409,
  187933,
  227047,
  250167,
  430745,
  285839,
  44841,
  215175,
  181493,
  234721,
  467439,
  457019,
  509077,
  359201,
  11027,
  429903,
  268385,
  97401,
  139475,
  412901,
  384039,
  216497,
  22511,
  368699,
  428213,
  124377,
  436615,
  360837,
  429077,
  400709,
  514399,
  241011,
  80915,
  66531,
  173107,
  524193,
  65449,
  357243,
  472359,
  457421,
  276537,
  332615,
  121161,
  30523,
  66625,
  493851,
  75357,
  72053,
  138243,
  498569,
  241009,
  482601,
  166199,
  302657,
  272589,
  92643,
  492011,
  118737,
  308979,
  6617,
  305097,
  466585,
  318419,
  149945,
  380153,
  291261,
  259013,
  88595,
  476949,
  41391,
  170045,
  256319,
  295207,
  372873,
  438169,
  196221,
  49225,
  200173,
  269399,
  202607,
  299353,
  221825,
  391325,
  248509,
  183861,
  62937,
  370937,
  308047,
  208963,
  362413,
  230899,
  518307,
  187739,
  460067,
  391577,
  434663,
  35647,
  286715,
  294589,
  389463,
  162933,
  411537,
  167093,
  65127,
  121753,
  344463,
  516765,
  414291,
  128981,
  243439,
  470701,
  437981,
  342493,
  322629,
  504431,
  383423,
  468615,
  402015,
  262491,
  225053,
  195803,
  145011,
  460657,
  289933,
  490493,
  258547,
  51845,
  284765,
  300867,
  142667,
  50691,
  172085,
  353993,
  500573,
  254085,
  113615,
  519945,
  129487,
  277663,
  164767,
  416047,
  155579,
  248353,
  109161,
  382985,
  388793,
  507679,
  299015,
  274125,
  204053,
  431273,
  513631,
  483277,
  36177,
  471323,
  306393,
  17351,
  87967,
  382551,
  4019,
  377623,
  114051,
  355573,
  523977,
  484437,
  271457,
  199003,
  108197,
  271007,
  344501,
  234951,
  157225,
  409627,
  129345,
  96097,
  264835,
  255921,
  338915,
  251267,
  277093,
  73221,
  392169,
  183719,
  139877,
  344877,
  195725,
  163315,
  6783,
  43537,
  244069,
  362533,
  97361,
  189961,
  240165,
  118595,
  452911,
  427757,
  314895,
  494405,
  21653,
  399265,
  330607,
  369479,
  491155,
  31265,
  434459,
  288753,
  84455,
  293565,
  460293,
  21425,
  61805,
  437609,
  259113,
  305619,
  234341,
  408429,
  294183,
  424917,
  521281,
  437047,
  430353,
  81257,
  393247,
  382355,
  131613,
  299665,
  20109,
  201617,
  185989,
  428861,
  94113,
  511567,
  412233,
  429035,
  519411,
  314931,
  519951,
  77705,
  87461,
  333213,
  284637,
  7935,
  371921,
  237783,
  428329,
  442883,
  322605,
  277795,
  103027,
  85391,
  176447,
  43693,
  323725,
  311337,
  344439,
  208419,
  93071,
  25893,
  395859,
  451291,
  448075,
  173979,
  397339,
  491231,
  93091,
  454985,
  351989,
  302041,
  361357,
  67663,
  170813,
  387243,
  224899,
  20489,
  484965,
  313779,
  361233,
  469957,
  332483,
  62329,
  246513,
  168613,
  11503,
  177541,
  333063,
  195349,
  295523,
  259807,
  142231,
  269037,
  482291,
  299483,
  252479,
  283629,
  240193,
  320163,
  486023,
  520959,
  63111,
  10213,
  431655,
  439659,
  147217,
  411661,
  308157,
  321907,
  501265,
  419839,
  371177,
  513683,
  419235,
  162683,
  346933,
  37079,
  455955,
  491201,
  343417,
  111105,
  516529,
  47713,
  171791,
  388113,
  247363,
  484029,
  169773,
  250871,
  439705,
  249861,
  93563,
  276649,
  347735,
  275159,
  444463,
  444457,
  340165,
  238219,
  511421,
  3711,
  293153,
  433595,
  498467,
  328597,
  462993,
  242767,
  280633,
  188319,
  147243,
  59519,
  10371,
  109385,
  421289,
  519729,
  311493,
  55557,
  62243,
  332767,
  132337,
  502769,
  190475,
  220289,
  368661,
  101033,
  159287,
  346783,
  88817,
  322189,
  297843,
  22851,
  314555,
  30719,
  226467,
  316875,
  64935,
  348957,
  411789,
  321345,
  100427,
  196969,
  216519,
  212671,
  238033,
  470295,
  428273,
  85607,
  368793,
  146055,
  169169,
  38043,
  26807,
  146255,
  400421,
  174601,
  491897,
  233119,
  387201,
  348875,
  68805,
  187145,
  20447,
  486529,
  60531,
  124855,
  456157,
  471813,
  183515,
  375031,
  107323,
  448423,
  283509,
  379319,
  417789,
  209295,
  238877,
  395877,
  427147,
  23411,
  223803,
  22595,
  108801,
  322185,
  58469,
  483411,
  394163,
  413329,
  354901,
  522441,
  109481,
  168295,
  269339,
  238611,
  503677,
  437089,
  379373,
  57093,
  91041,
  125787,
  324727,
  177279,
  50247,
  416961,
  506185,
  404993,
  374955,
  52633,
  301157,
  236643,
  167913,
  393529,
  178293,
  73215,
  385523,
  260125,
  486055,
  449495,
  58677,
  375277,
  320979,
  451277,
  330395,
  438251,
  310281,
  317481,
  126541,
  214215,
  475315,
  315097,
  318983,
  317505,
  474263,
  496911,
  125405,
  43067,
  33517,
  389253,
  42313,
  129447,
  421007,
  362237,
  227073,
  126017,
  479819,
  210921,
  79791,
  198363,
  384549,
  329489,
  44665,
  356461,
  134445,
  135127,
  449099,
  155947,
  168329,
  365601,
  353765,
  419155,
  150437,
  119517,
  298069,
  28947,
  35443,
  413727,
  269677,
  366217,
  136227,
  155333,
  464397,
  130865,
  241509,
  320225,
  137493,
  271105,
  75223,
  375527,
  143779,
  112229,
  245705,
  112867,
  411037,
  89325,
  237307,
  161365,
  484761,
  207255,
  206127,
  321743,
  384985,
  442877,
  263129,
  211505,
  325823,
  85081,
  445183,
  409583,
  283329,
  65059,
  131307,
  243747,
  387529,
  423355,
  2995,
  491469,
  125507,
  411969,
  27743,
  437423,
  506305,
  259567,
  3309,
  286473,
  358905,
  160297,
  4615,
  42525,
  152491,
  237177,
  171253,
  356853,
  320371,
  66075,
  134765,
  365639,
  276209,
  498613,
  353187,
  176017,
  29435,
  113007,
  501773,
  207885,
  370521,
  345009,
  148935,
  43737,
  36005,
  386465,
  137591,
  20111,
  52163,
  84625,
  67001,
  244925,
  408647,
  90653,
  112981,
  416825,
  192987,
  291471,
  408537,
  370213,
  487317,
  195075,
  484829,
  452531,
  491715,
  175131,
  296705,
  182743,
  401533,
  145133,
  47655,
  106901,
  385389,
  333105,
  495649,
  198013,
  490969,
  15769,
  517559,
  360277,
  99815,
  523297,
  61605,
  140317,
  140553,
  206363,
  373783,
  441221,
  331555,
  53047,
  15759,
  252725,
  100597,
  410343,
  471953,
  159549,
  327959,
  441471,
  110655,
  283353,
  356177,
  63757,
  264055,
  154289,
  234919,
  175549,
  66581,
  218069,
  47281,
  185375,
  94719,
  117301,
  486893,
  239833,
  54185,
  320917,
  67735,
  104983,
  304465,
  440313,
  470631,
  28787,
  381879,
  207727,
  98081,
  119355,
  141793,
  522435,
  499393,
  95157,
  251465,
  93685,
  247587,
  462773,
  360895,
  193175,
  47033,
  421069,
  240597,
  9847,
  496635,
  353401,
  21525,
  514603,
  215633,
  324941,
  266803,
  195843,
  457817,
  297245,
  87983,
  247833,
  404853,
  209401,
  130941,
  371029,
  274429,
  141599,
  340233,
  402543,
  243355,
  502569,
  219343,
  316665,
  409117,
  467413,
  167563,
  8767,
  490011,
  31969,
  244789,
  173937,
  488819,
  279215,
  331,
  31395,
  79851,
  360451,
  91743,
  453003,
  247647,
  50239,
  130511,
  266153,
  436609,
  305001,
  386119,
  51055,
  391537,
  118293,
  189295,
  436559,
  6587,
  484879,
  185985,
  369453,
  39209,
  85417,
  61079,
  515227,
  378617,
  369691,
  289961,
  417129,
  216939,
  66181,
  161249,
  221775,
  523799,
  291817,
  170383,
  418119,
  363171,
  454299,
  207721,
  447149,
  456603,
  391507,
  36489,
  152227,
  75227,
  172843,
  68107,
  452955,
  145421,
  150221,
  246009,
  78407,
  287993,
  399961,
  344161,
  224507,
  296827,
  435787,
  373389,
  405083,
  299121,
  224369,
  86133,
  289307,
  341689,
  380379,
  247229,
  111781,
  193173,
  381567,
  90409,
  87017,
  379085,
  195549,
  328503,
  448819,
  521439,
  523211,
  317559,
  508697,
  435439,
  31307,
  300425,
  516943,
  337803,
  476487,
  78833,
  168493,
  82991,
  1813,
  279039,
  139735,
  212037,
  507637,
  27191,
  342299,
  333581,
  217225,
  361105,
  141349,
  217139,
  455297,
  22027,
  208285,
  105195,
  343399,
  365187,
  166003,
  267747,
  173265,
  104309,
  263875,
  437963,
  390657,
  186295,
  221165,
  285495,
  385165,
  350075,
  373929,
  47901,
  362017,
  223307,
  221823,
  171781,
  317267,
  452633,
  120487,
  425947,
  71937,
  144669,
  30537,
  320807,
  496709,
  36143,
  35249,
  54399,
  457465,
  374257,
  259027,
  434899,
  192051,
  107773,
  215557,
  272847,
  145053,
  179531,
  149193,
  455445,
  7111,
  146859,
  325841,
  500423,
  330427,
  332067,
  224651,
  142941,
  281513,
  524213,
  43221,
  44235,
  440259,
  330081,
  413537,
  455371,
  58199,
  47729,
  95937,
  211567,
  206227,
  427361,
  466909,
  1005,
  206647,
  104791,
  257555,
  255555,
  236197,
  24417,
  47507,
  357477,
  347973,
  131241,
  452411,
  494475,
  342379,
  489111,
  6675,
  509787,
  329523,
  142947,
  212759,
  502533,
  12323,
  499179,
  353173,
  5519,
  146257,
  136979,
  333593,
  242981,
  76435,
  369869,
  520221,
  170155,
  299129,
  373969,
  297759,
  168945,
  258359,
  383557,
  73439,
  304311,
  488225,
  185907,
  242127,
  479013,
  296237,
  143831,
  376901,
  384261,
  70843,
  267839,
  390079,
  207373,
  336845,
  188153,
  441483,
  502729,
  144489,
  143595,
  521973,
  347585,
  174495,
  362967,
  195351,
  430933,
  331571,
  155305,
  78807,
  252037,
  67461,
  105441,
  204725,
  517849,
  236481,
  441255,
  272255,
  396343,
  120681,
  105631,
  236671,
  192703,
  130503,
  470853,
  156907,
  370583,
  346321,
  412653,
  332919,
  284217,
  395281,
  283305,
  108461,
  518885,
  321629,
  465913,
  138519,
  499083,
  391627,
  318493,
  146001,
  167501,
  478239,
  469563,
  311423,
  82501,
  337647,
  462381,
  496957,
  417031,
  449591,
  322099,
  110997,
  361467,
  124215,
  310941,
  457079,
  343763,
  274689,
  454023,
  196851,
  140011,
  329071,
  387481,
  408773,
  40235,
  351537,
  463559,
  73109,
  366319,
  370429,
  61679,
  399063,
  252481,
  429989,
  58691,
  470643,
  485623,
  402333,
  225503,
  144559,
  121153,
  76475,
  350553,
  3935,
  31375,
  94239,
  302415,
  401575,
  211757,
  271075,
  159559,
  463171,
  251429,
  357403,
  432513,
  106507,
  255729,
  429641,
  291489,
  71215,
  99219,
  45203,
  83541,
  203975,
  337455,
  142077,
  129255,
  466479,
  110353,
  222145,
  184543,
  4223,
  61693,
  480821,
  124743,
  487855,
  176403,
  194931,
  380825,
  487637,
  337069,
  508087,
  497721,
  181339,
  398547,
  516409,
  139969,
  185615,
  25049,
  51321,
  279923,
  195797,
  80977,
  520641,
  206053,
  51091,
  226235,
  475763,
  449191,
  25233,
  390933,
  365715,
  114355,
  336241,
  345473,
  253533,
  170419,
  276927,
  52391,
  192871,
  348987,
  404949,
  365431,
  150615,
  363203,
  273999,
  25469,
  341535,
  197597,
  491999,
  191295,
  360335,
  436757,
  400821,
  488647,
  292951,
  164863,
  419047,
  410339,
  364887,
  447227,
  312255,
  401577,
  491153,
  520151,
  28657,
  137485,
  288963,
  172713,
  366519,
  293067,
  505061,
  221561,
  376689,
  209791,
  274969,
  61145,
  143303,
  157459,
  518153,
  349473,
  169179,
  3579,
  364721,
  142975,
  311275,
  474695,
  85031,
  14031,
  57731,
  399487,
  51707,
  103261,
  38863,
  356271,
  472975,
  160727,
  168239,
  468851,
  130069,
  188281,
  515097,
  119167,
  74835,
  62219,
  81827,
  573,
  464349,
  323045,
  201719,
  332681,
  131117,
  288875,
  509111,
  341827,
  59475,
  352861,
  83495,
  3797,
  30739,
  454485,
  411881,
  167767,
  168457,
  268881,
  425747,
  186685,
  351955,
  108987,
  475831,
  222691,
  399295,
  422525,
  353383,
  397973,
  200631,
  91443,
  267695,
  217305,
  454703,
  479915,
  357239,
  383125,
  425711,
  449165,
  505985,
  470431,
  12239,
  514817,
  348329,
  304277,
  465579,
  128661,
  238205,
  432041,
  30891,
  324155,
  486707,
  189683,
  124635,
  76543,
  335841,
  330153,
  42873,
  33289,
  404423,
  432881,
  372293,
  56239,
  22533,
  479905,
  378149,
  245273,
  141797,
  343875,
  146241,
  271381,
  297861,
  424403,
  36713,
  211755,
  387069,
  243775,
  216099,
  258397,
  100383,
  459721,
  415609,
  330499,
  417793,
  87947,
  314485,
  5171,
  472093,
  133175,
  517795,
  57461,
  159591,
  518095,
  385625,
  477059,
  301889,
  278983,
  449397,
  430415,
  484097,
  214421,
  64161,
  484031,
  252913,
  147907,
  497375,
  79797,
  321555,
  36029,
  403759,
  306137,
  106975,
  278317,
  352455,
  311517,
  257131,
  98943,
  513449,
  68961,
  295611,
  9281,
  106999,
  37679,
  410957,
  356071,
  237009,
  452951,
  120029,
  4327,
  299495,
  302453,
  266965,
  75049,
  299435,
  33105,
  105635,
  58447,
  468661,
  49143,
  184387,
  304455,
  161973,
  383407,
  161591,
  99557,
  136555,
  227539,
  53703,
  183421,
  1309,
  219647,
  181273,
  225657,
  499885,
  7793,
  316853,
  203725,
  510555,
  439029,
  140139,
  226933,
  53193,
  426793,
  323911,
  197289,
  451239,
  400865,
  186897,
  415817,
  247495,
  338897,
  81519,
  508191,
  51431,
  138341,
  168981,
  509307,
  425997,
  31671,
  466107,
  107807,
  119149,
  45445,
  463293,
  169005,
  494957,
  233889,
  379071,
  56207,
  370557,
  190641,
  313037,
  409251,
  115275,
  432787,
  149865,
  128547,
  387723,
  208827,
  513141,
  342985,
  451773,
  123373,
  266837,
  356005,
  210101,
  222037,
  309733,
  291245,
  198909,
  109053,
  338203,
  358821,
  315763,
  260497,
  270577,
  50711,
  16659,
  432321,
  426009,
  363977,
  169873,
  397129,
  506275,
  34723,
  406577,
  512633,
  446285,
  126645,
  359527,
  234217,
  248385,
  72351,
  490133,
  353629,
  275083,
  257899,
  286679,
  506675,
  392463,
  64919,
  304447,
  167823,
  469789,
  128867,
  341797,
  269127,
  139817,
  73833,
  474985,
  146207,
  259459,
  514783,
  225979,
  379661,
  65709,
  159455,
  68921,
  406531,
  165331,
  162783,
  134373,
  524181,
  325319,
  382207,
  374607,
  459885,
  110169,
  109615,
  207343,
  242825,
  26893,
  304861,
  271033,
  39067,
  69399,
  398623,
  82457,
  68137,
  30655,
  357523,
  110607,
  70321,
  200023,
  333595,
  519989,
  406867,
  454885,
  371469,
  181747,
  287901,
  328919,
  452885,
  268871,
  114371,
  174289,
  204689,
  339681,
  244961,
  138059,
  23775,
  157431,
  364049,
  168241,
  386843,
  304593,
  244153,
  455883,
  135783,
  189137,
  72541,
  436449,
  17323,
  184745,
  281281,
  392061,
  349459,
  451941,
  36721,
  442071,
  215979,
  96757,
  148521,
  45921,
  77641,
  362129,
  362027,
  1419,
  213361,
  206533,
  416007,
  304031,
  432157,
  55333,
  318205,
  399399,
  259627,
  269895,
  113389,
  470083,
  267011,
  224557,
  4671,
  387213,
  73547,
  288035,
  267849,
  439837,
  294697,
  236133,
  59493,
  188573,
  453863,
  27265,
  422439,
  183715,
  310065,
  436519,
  347411,
  378187,
  92151,
  444977,
  215793,
  388667,
  346607,
  219657,
  470907,
  105075,
  382637,
  318467,
  193169,
  405245,
  304635,
  140987,
  100999,
  44247,
  193987,
  333067,
  504099,
  458487,
  78945,
  113979,
  327675,
  182731,
  190221,
  109961,
  47721,
  137905,
  509935,
  415927,
  296613,
  151271,
  397933,
  475833,
  209163,
  122023,
  121877,
  7681,
  390363,
  188591,
  37265,
  102951,
  510389,
  491125,
  184595,
  150425,
  174379,
  496253,
  397359,
  142605,
  32063,
  278813,
  176873,
  85853,
  234275,
  234797,
  377841,
  484883,
  342529,
  218631,
  200297,
  328805,
  400701,
  223069,
  181567,
  309437,
  405825,
  205209,
  296555,
  308911,
  411577,
  410949,
  13167,
  272245,
  163751,
  136507,
  110507,
  165829,
  435515,
  63761,
  474227,
  445641,
  163547,
  311013,
  18131,
  189441,
  78087,
  390825,
  491063,
  154625,
  445825,
  178395,
  160541,
  129953,
  17663,
  59561,
  259215,
  418625,
  345359,
  268621,
  225865,
  102075,
  372677,
  227767,
  437373,
  403481,
  151031,
  477203,
  150857,
  147693,
  245381,
  231701,
  201019,
  284969,
  41339,
  430169,
  7919,
  243345,
  9655,
  391383,
  345769,
  58521,
  329835,
  178569,
  474807,
  66627,
  231369,
  439087,
  298597,
  200877,
  512175,
  352351,
  87387,
  269711,
  458703,
  369055,
  402877,
  252607,
  197639,
  213247,
  232615,
  92657,
  429865,
  91295,
  60675,
  515225,
  394623,
  167027,
  91441,
  168589,
  407489,
  61537,
  12727,
  419841,
  116529,
  431501,
  27431,
  59541,
  422421,
  352663,
  160503,
  105091,
  440027,
  391909,
  499147,
  227823,
  356185,
  209859,
  237757,
  294719,
  198541,
  93801,
  193925,
  226339,
  449107,
  131637,
  385919,
  191353,
  507939,
  119171,
  120225,
  194801,
  420895,
  390275,
  357333,
  82391,
  369141,
  57705,
  477633,
  468597,
  497491,
  329805,
  363165,
  314103,
  41947,
  217175,
  456325,
  508435,
  370641,
  149311,
  377971,
  492869,
  88373,
  220541,
  269265,
  519127,
  302535,
  450231,
  114681,
  439607,
  48657,
  195225,
  441259,
  270287,
  123281,
  157339,
  24067,
  277467,
  23209,
  480037,
  75733,
  120049,
  74375,
  308017,
  289173,
  95853,
  45987,
  115339,
  466973,
  506551,
  254681,
  334965,
  243101,
  238815,
  326495,
  172489,
  73113,
  119345,
  79275,
  119931,
  30161,
  240745,
  428825,
  47675,
  121967,
  233195,
  123165,
  479149,
  97037,
  142089,
  257253,
  214405,
  251141,
  119029,
  514015,
  248751,
  517863,
  142393,
  154935,
  12367,
  106931,
  267209,
  416673,
  176893,
  384475,
  119825,
  164611,
  358501,
  357389,
  452533,
  460727,
  142831,
  31497,
  408089,
  392383,
  67637,
  288857,
  175585,
  266747,
  41765,
  393901,
  72789,
  107687,
  30391,
  285671,
  213485,
  146743,
  428499,
  168767,
  73743,
  253675,
  333159,
  40307,
  454903,
  474825,
  437973,
  44739,
  308231,
  179883,
  134153,
  475069,
  188537,
  508973,
  220637,
  188005,
  192269,
  383763,
  408723,
  173275,
  34577,
  399119,
  84701,
  198817,
  328407,
  198645,
  90179,
  468431,
  116563,
  333127,
  211005,
  246461,
  268813,
  96113,
  3447,
  212257,
  39265,
  168843,
  281933,
  149947,
  392735,
  194069,
  418871,
  378107,
  399309,
  45741,
  396531,
  448939,
  388199,
  474421,
  41241,
  213887,
  56587,
  177657,
  494225,
  278607,
  176267,
  59741,
  270589,
  165791,
  155361,
  6473,
  162547,
  115111,
  453061,
  50925,
  405033,
  489907,
  471277,
  514717,
  267253,
  417221,
  300903,
  482243,
  48449,
  322747,
  130635,
  506133,
  206951,
  93809,
  511615,
  106977,
  511393,
  165785,
  289257,
  105363,
  293217,
  30345,
  79175,
  60301,
  461261,
  401323,
  5943,
  491845,
  148305,
  104025,
  39305,
  474573,
  341611,
  274603,
  441219,
  334337,
  195011,
  183303,
  96023,
  495489,
  380921,
  203041,
  388783,
  235133,
  281039,
  38607,
  62285,
  123401,
  89741,
  410791,
  32253,
  173905,
  342407,
  315951,
  33885,
  336275,
  243533,
  283835,
  290881,
  183377,
  305311,
  467285,
  464893,
  197217,
  61965,
  45613,
  517645,
  313339,
  450115,
  167999,
  514009,
  68119,
  32235,
  77667,
  29657,
  462767,
  231101,
  236467,
  173483,
  378741,
  437545,
  428441,
  487347,
  125753,
  439559,
  313519,
  176115,
  203733,
  361193,
  66291,
  372245,
  228819,
  457165,
  353889,
  271861,
  2961,
  290687,
  114523,
  106185,
  115973,
  301141,
  143421,
  471109,
  499515,
  398857,
  512545,
  331319,
  159897,
  309705,
  426123,
  196201,
  253965,
  441231,
  363183,
  413997,
  72719,
  418185,
  136949,
  84053,
  66205,
  385911,
  166495,
  440155,
  428225,
  19885,
  56283,
  138471,
  26285,
  381021,
  465113,
  266033,
  496847,
  419331,
  494525,
  10463,
  460669,
  63139,
  441583,
  144601,
  151503,
  503305,
  84051,
  24255,
  497473,
  97709,
  184819,
  263255,
  115365,
  308643,
  199949,
  261511,
  227647,
  83109,
  144957,
  254307,
  102661,
  95753,
  392433,
  310997,
  203229,
  308719,
  240017,
  138373,
  29295,
  254391,
  16235,
  143789,
  257247,
  265947,
  92127,
  110075,
  164045,
  455097,
  122351,
  305117,
  137161,
  84205,
  512991,
  422887,
  442063,
  75573,
  292001,
  141359,
  14853,
  306481,
  102725,
  81015,
  14697,
  86609,
  146273,
  238385,
  385869,
  485749,
  60087,
  147359,
  62247,
  117087,
  352319,
  38317,
  404431,
  361969,
  148009
};

} // namespace Dakota

#endif