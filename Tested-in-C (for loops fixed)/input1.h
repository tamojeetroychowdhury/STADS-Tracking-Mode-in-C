#include <stdio.h>
#include <stdlib.h>

;int len1 = 26;
int len2 = 24;
int len3;

double p1[40][3];
double p2[40][3];
int comm;
double pred[40][3];
double matchmat[40][3];

double truemat[40][2];

int good = 0;
int num = 18;
double r = num*pixel_size;

double e = 0.000001;
double snt_out[5][3], snt_match=0;



double prev1[40][3] =
{
-6.27771192990944e-05 , 6.28869545683172e-05 , 31 ,
0.000650077897040617 , -0.000403966053429015 , 32 ,
-0.000722298931579938 , 0.00100747125328688 , 70 ,
-0.00190817408476628 , -0.000687001448112021 , 252 ,
0.000333208385475908 , -0.000816512890427677 , 409 ,
0.000650192437312053 , -0.000404372584450126 , 633 ,
-0.00108076129700211 , 0.000133309459601561 , 1129 ,
0.00066478319681518 , 0.000107420482825636 , 1500 ,
-0.00190620969699828 , 0.00138286450213734 , 1745 ,
-0.000616658537898479 , -0.000181261209292924 , 2377 ,
-0.000485128318345485 , 9.23820326606437e-05 , 2436 ,
-0.00113551872559765 , -0.00134714878659792 , 4066 ,
-0.000400085878209512 , -0.000105166923051751 , 4530 ,
0.000705801540900375 , 0.00102924596638064 , 4655 ,
0.000453302493217608 , -0.00142353567066423 , 4915 ,
-0.00170124286120857 , 0.0011486310556729 , 5948 ,
-0.000624090293466841 , 0.000826722209435255 , 6136 ,
0.000627322841854635 , -0.000958801126836505 , 6310 ,
-0.000400077806196071 , 0.000167873091029623 , 6542 ,
-0.000332561186867334 , 0.000826460436333273 , 6726 ,
0.00102704048613991 , -0.00019791618520349 , 7059 ,
-0.00120551717417077 , -0.00125967577853486 , 7776 ,
-0.000161553215838179 , -0.00122606311809489 , 7955 ,
0.000792185630781686 , -0.0010040496842578 , 8040 ,
-0.00172973655849645 , -0.000117039585116407 , 8122 ,
-0.00108563724087431 , 0.00132496113668565 , 8804
};


double prev2[40][3] =
{
-0.000141174088116146 , 0.000141590934807836 , 31 ,
0.00057062617610829 , -0.000326795363799492 , 32 ,
-0.000798716072300519 , 0.0010877250600902 , 70 ,
0.000252874750460045 , -0.000738621299795706 , 409 ,
0.00057073980346805 , -0.000327202108052889 , 633 ,
-0.00115908884497346 , 0.000214246031857276 , 1129 ,
0.000586464673035439 , 0.000184531999526818 , 1500 ,
-0.000695605940520748 , -0.000101354252584696 , 2377 ,
-0.000563479559238399 , 0.000172011139513134 , 2436 ,
-0.00121697812367594 , -0.0012660928484502 , 4066 ,
-0.000478857298215181 , -2.57295895872426e-05 , 4530 ,
0.000629527861777997 , 0.0011062958309554 , 4655 ,
0.000371618825041509 , -0.00134580422565187 , 4915 ,
-0.00177757522818189 , 0.00123111641735157 , 5948 ,
-0.00070087999567384 , 0.000906727757069652 , 6136 ,
0.000546643602985499 , -0.000881517042438554 , 6310 ,
-0.000478260942491579 , 0.000247319390375763 , 6542 ,
-0.000409318577505941 , 0.000905814698596577 , 6726 ,
0.00094800465971362 , -0.000121578630555832 , 7059 ,
-0.00128679698227178 , -0.00117847967734121 , 7776 ,
-0.000242747857742776 , -0.00114707120293969 , 7955 ,
0.000711382693288981 , -0.000927108585701155 , 8040 ,
-0.00180870249992042 , -3.47082823853803e-05 , 8122 ,
-0.00116144710366268 , 0.00140609794699343 , 8804
};

int arr_len3[] = {23, 22, 21, 22, 22, 23, 21};
double arr_truemat[200][2] =
{
-0.000219482425261588 , 0.000220407481186828 ,
0.000491272430404318 , -0.00024951938769284 ,
-0.000875055384266258 , 0.00116810354475697 ,
0.000172638432073634 , -0.000660627170432887 ,
0.000491385148286462 , -0.00024992634941659 ,
-0.00123733802744653 , 0.000295300550931328 ,
0.000508241377901931 , 0.000261753487562048 ,
-0.000774468833388462 , -2.13346514947489e-05 ,
-0.000641746384042191 , 0.00025175497232365 ,
-0.00129835216606036 , -0.00118493387454311 ,
-0.00055754247210707 , 5.38201670903275e-05 ,
0.000553344923406384 , 0.00118346379162247 ,
0.000290037326880904 , -0.00126797660155736 ,
-0.00185383978748222 , 0.00131373316079755 ,
-0.000777590128138779 , 0.000986855760928028 ,
0.000466065402714505 , -0.000804132853671228 ,
-0.000556359217928505 , 0.000326880758549083 ,
-0.00048599365877094 , 0.000985289982689134 ,
0.000869069146407355 , -4.51352705979715e-05 ,
-0.00136799265030337 , -0.00109717932037051 ,
-0.000323847727554021 , -0.00106797892733254 ,
0.00063068267692201 , -0.000850068385821066 ,
-0.00188759548902936 , 4.77416736568127e-05 ,




-0.000297703394275344 , 0.000299337745780745 ,
0.000412015366065286 , -0.000172137041423887 ,
-0.000951318060276175 , 0.00124860792747936 ,
9.24980894830584e-05 , -0.000582529384986637 ,
0.000412127177865934 , -0.000172544224864657 ,
-0.00131551012734044 , 0.000376474273665683 ,
0.000430112067696377 , 0.000339086026617773 ,
-0.000853248516608152 , 5.87988030340464e-05 ,
-0.000719930063658219 , 0.000331614726046254 ,
-0.0013796422789796 , -0.00110367059684318 ,
-0.000636142687394014 , 0.000133483533689949 ,
0.000477251569973768 , 0.00126075092226269 ,
0.000208556598592667 , -0.00119005168706585 ,
-0.00193003774912617 , 0.00139648261280669 ,
-0.000854221897355914 , 0.00106710743043601 ,
0.000385586891350737 , -0.000726647470503317 ,
-0.000634373894178957 , 0.000406558381661054 ,
-0.000562587629493727 , 0.00106488746766556 ,
0.000790232679525739 , 3.14149407436851e-05 ,
-0.00144910559767337 , -0.00101577343303429 ,
-0.000404854217025664 , -0.000988785122040098 ,
0.000550084230240795 , -0.000772928009980664 ,




-0.000375838254740112 , 0.000378382887909604 ,
0.000332853693598955 , -9.46472345012972e-05 ,
-0.00102750528971401 , 0.00132923943630111 ,
1.2452386316335e-05 , -0.000504326819288439 ,
0.000332964602676497 , -9.50546439145118e-05 ,
-0.00139360642399018 , 0.000457768464683581 ,
0.000352075502855309 , 0.000416530703314166 ,
-0.000931946286417953 , 0.000139047327496962 ,
-0.000798031865179375 , 0.000411591603092415 ,
-0.00146084988450217 , -0.00102230174003253 ,
-0.000714659227671621 , 0.000213261704276217 ,
0.000401246649479891 , 0.00133815830389064 ,
0.000127175244833552 , -0.0011120283642804 ,
-0.00093077650631497 , 0.00114748398276496 ,
0.000305206723827282 , -0.000649059796244747 ,
-0.000712306229106714 , 0.000486353453260796 ,
-0.000639101684909935 , 0.00114460834017103 ,
0.00071149399675928 , 0.000108073056254667 ,
-0.00153013723969805 , -0.000934260733392154 ,
-0.00048576871342666 , -0.0009094886109395 ,
0.000469586006503403 , -0.000695686376958722 ,




-0.000453888262341613 , 0.000457544074201856 ,
0.000253786127822629 , -1.70488695702423e-05 ,
-0.00110361825860363 , 0.00140999930718789 ,
-6.75000094031912e-05 , -0.000426018342293331 ,
0.000253896137497944 , -1.74565092202417e-05 ,
-0.0014716281932103 , 0.000539184396463868 ,
0.000274130447912734 , 0.000494088611262097 ,
-0.00101056343525197 , 0.000219412145896926 ,
-0.000876053051951689 , 0.00049168681339909 ,
-0.00154197640061309 , -0.000940826021433249 ,
-0.000793093372685176 , 0.000293155880336615 ,
0.000325329013685974 , 0.00141568702471097 ,
4.58918749887804e-05 , -0.00103390550866207 ,
-0.00100725515455752 , 0.00122798664290137 ,
0.000224923559634663 , -0.000571368727488686 ,
-0.000790157476827892 , 0.000566267174402271 ,
-0.000715537016714685 , 0.00122445379450764 ,
0.00063285184012572 , 0.000184840135486275 ,
-0.00161108898767828 , -0.000852639932075465 ,
0.00146325613802353 , -0.00138787377953849 ,
-0.000566592599612932 , -0.000830088210952855 ,
0.000389186663599532 , -0.000618342398893226 ,




-0.000531854668928445 , 0.000536822478659031 ,
0.000174811387801959 , 6.06591576456439e-05 ,
-0.00014736042532269 , -0.000347602816018663 ,
0.00017492050135862 , 6.02512834854483e-05 ,
-0.00154957670735819 , 0.000620723349416174 ,
0.000196275671443208 , 0.000571760851120811 ,
-0.00108910125180426 , 0.000299894489814853 ,
-0.000953994883635548 , 0.000571901574497615 ,
-0.00162302324129016 , -0.000859242150931008 ,
-0.000871446398392138 , 0.000373167270849591 ,
-3.52948968989021e-05 , -0.000955681988969366 ,
-0.00108365903823298 , 0.00130861664371708 ,
0.000144736062754743 , -0.000493573154052957 ,
-0.000867928887771225 , 0.000646300753712287 ,
-0.000791894813117474 , 0.00130442503270168 ,
0.000554304955911809 , 0.000261717244810115 ,
0.000981529533848575 , -0.00141495849983279 ,
-0.00169196224897532 , -0.000770909732210305 ,
0.00138105639488789 , -0.00131256924868169 ,
-0.000647327254098672 , -0.000750582731980917 ,
0.00030888486398455 , -0.000540894981224253 ,
-0.000670154616544446 , -0.00144931820871564 ,




-0.000609738722571973 , 0.000616219282721142 ,
9.5928196789742e-05 , 0.000138477958403016 ,
-0.000227130184825685 , -0.000269079095482133 ,
9.60364174742271e-05 , 0.000138069845449983 ,
-0.00162745323539748 , 0.000702386611956963 ,
0.000118509946004897 , 0.000649548530656322 ,
-0.00116756102109283 , 0.000380495598480784 ,
-0.00103185861626696 , 0.000652237111583389 ,
-0.00170399181658026 , -0.000777548830897809 ,
-0.000949719577024676 , 0.000453297092353428 ,
-0.0001163864522001 , -0.000877356667196935 ,
-0.00115998935015474 , 0.00138937522604034 ,
6.46429015946528e-05 , -0.00041567195892073 ,
-0.000945621708738393 , 0.00072645540745928 ,
-0.000775826612458758 , -0.00140716365714308 ,
-0.000868176258897126 , 0.00138452326457154 ,
0.000475852094615082 , 0.000338705457473305 ,
0.000899299811563238 , -0.00133869915444376 ,
-0.00177275842708662 , -0.000689068829339502 ,
0.00129896814326048 , -0.00123717107324792 ,
-0.000727974051127623 , -0.00067097097683603 ,
0.000228679274613576 , -0.000463343022636152 ,
-0.000752293827314916 , -0.00136959718594213 ,




-0.000687541667625862 , 0.00069573567533286 ,
1.71352821651151e-05 , 0.000216408650998282 ,
-0.000306810607097482 , -0.000190446028639086 ,
1.72426131870108e-05 , 0.000216000294960379 ,
-0.00170525904296106 , 0.00078417548058644 ,
4.08320480829856e-05 , 0.000727452764800514 ,
-0.00124594402452333 , 0.00046121671884584 ,
-0.00110964550231849 , 0.000732694657586379 ,
-0.00178488353267518 , -0.000695744756113897 ,
-0.00102791417715185 , 0.000533546569015911 ,
-0.000197384167769219 , -0.000798928398513574 ,
-1.53572510788638e-05 , -0.000337664018180515 ,
-0.00102323718296401 , 0.000806732359622917 ,
-0.00085795600097707 , -0.00132728856383739 ,
0.000397492010886005 , 0.000415805853654375 ,
0.000817176916768941 , -0.00126234490425229 ,
-0.00185347892172119 , -0.000607115911343838 ,
0.00121698999317072 , -0.00116167823981453 ,
-0.000808534360743945 , -0.000591251741174296 ,
0.000148568566875983 , -0.000385685414999036 ,
-0.000834342912094459 , -0.00128977571691641
}


;

