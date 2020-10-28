#!/usr/bin/env python

import os
import sys
'''samples output:
	#local atom-atom
   Warning==> Structures have not yet been fitted.
   Fitting structures...
   RMS: 0.722
   #local atom-atom by residue
    H95  ALA  :     H95  ALA    RMS: 0.267
    H96  GLU  :     H96  GLU    RMS: 0.612
    H97  ARG  :     H97  ARG    RMS: 1.024
    H98  LEU  :     H98  LEU    RMS: 0.612
    H99  ARG  :     H99  ARG    RMS: 0.558
   H100  ARG  :    H100  ARG    RMS: 1.135
   H100A THR  :    H100A THR    RMS: 0.489
   H100B PHE  :    H100B PHE    RMS: 0.338
   H101  ASP  :    H101  ASP    RMS: 0.952
   H102  TYR  :    H102  TYR    RMS: 0.456
   RMS: 1.727
   #local Ca-Ca
   Warning==> Structures have not yet been fitted.
   Fitting structures...
   RMS: 0.267
   #local Ca-Ca by residue
    H95  ALA  :     H95  ALA    RMS: 0.075
    H96  GLU  :     H96  GLU    RMS: 0.212
    H97  ARG  :     H97  ARG    RMS: 0.266
    H98  LEU  :     H98  LEU    RMS: 0.230
    H99  ARG  :     H99  ARG    RMS: 0.315
   H100  ARG  :    H100  ARG    RMS: 0.235
   H100A THR  :    H100A THR    RMS: 0.364
   H100B PHE  :    H100B PHE    RMS: 0.205
   H101  ASP  :    H101  ASP    RMS: 0.214
   H102  TYR  :    H102  TYR    RMS: 0.404
   RMS: 1.386
   #global atom-atom
   Warning==> Structures have not yet been fitted.
   Fitting structures...
   Warning==> Ignored reference atom  O   not found in mobile.
      Reference VAL    110  Mobile VAL    110
   Warning==> Ignored reference atom  O   not found in mobile.
      Reference SER    113  Mobile SER    113
   RMS: 1.551
   #global atom-atom by residue
     L1  GLU  :      L1  GLU    RMS: 2.861
     L2  LEU  :      L2  LEU    RMS: 1.577
     L3  VAL  :      L3  VAL    RMS: 0.930
     L4  MET  :      L4  MET    RMS: 0.881
     L5  THR  :      L5  THR    RMS: 0.415
     L6  GLN  :      L6  GLN    RMS: 0.320
     L7  THR  :      L7  THR    RMS: 2.057
     L8  PRO  :      L8  PRO    RMS: 1.901
     L9  LEU  :      L9  LEU    RMS: 1.458
    L10  SER  :     L10  SER    RMS: 0.876
    L11  LEU  :     L11  LEU    RMS: 1.316
    L12  PRO  :     L12  PRO    RMS: 0.463
    L13  VAL  :     L13  VAL    RMS: 1.145
    L14  SER  :     L14  SER    RMS: 1.075
    L15  LEU  :     L15  LEU    RMS: 1.763
    L16  GLY  :     L16  GLY    RMS: 1.847
    L17  ASP  :     L17  ASP    RMS: 1.491
    L18  GLN  :     L18  GLN    RMS: 1.657
    L19  ALA  :     L19  ALA    RMS: 0.868
    L20  SER  :     L20  SER    RMS: 0.944
    L21  ILE  :     L21  ILE    RMS: 0.651
    L22  SER  :     L22  SER    RMS: 1.036
    L23  CYS  :     L23  CYS    RMS: 0.516
    L24  ARG  :     L24  ARG    RMS: 4.115
    L25  SER  :     L25  SER    RMS: 0.490
    L26  SER  :     L26  SER    RMS: 1.019
    L27  GLN  :     L27  GLN    RMS: 0.668
    L28  SER  :     L28  SER    RMS: 0.459
    L29  LEU  :     L29  LEU    RMS: 0.459
    L30  VAL  :     L30  VAL    RMS: 0.410
    L30A HIS  :     L30A HIS    RMS: 1.247
    L30B SER  :     L30B SER    RMS: 0.760
    L30C ASN  :     L30C ASN    RMS: 1.047
    L30D GLY  :     L30D GLY    RMS: 0.831
    L30E ASN  :     L30E ASN    RMS: 1.178
    L31  THR  :     L31  THR    RMS: 0.353
    L32  TYR  :     L32  TYR    RMS: 1.471
    L33  LEU  :     L33  LEU    RMS: 0.389
    L34  HIS  :     L34  HIS    RMS: 1.141
    L35  TRP  :     L35  TRP    RMS: 0.392
    L36  TYR  :     L36  TYR    RMS: 1.444
    L37  LEU  :     L37  LEU    RMS: 0.443
    L38  GLN  :     L38  GLN    RMS: 1.110
    L39  LYS  :     L39  LYS    RMS: 1.237
    L40  PRO  :     L40  PRO    RMS: 1.359
    L41  GLY  :     L41  GLY    RMS: 1.221
    L42  GLN  :     L42  GLN    RMS: 1.224
    L43  SER  :     L43  SER    RMS: 0.656
    L44  PRO  :     L44  PRO    RMS: 0.454
    L45  LYS  :     L45  LYS    RMS: 0.823
    L46  LEU  :     L46  LEU    RMS: 0.877
    L47  LEU  :     L47  LEU    RMS: 1.476
    L48  ILE  :     L48  ILE    RMS: 0.272
    L49  TYR  :     L49  TYR    RMS: 1.478
    L50  LYS  :     L50  LYS    RMS: 0.978
    L51  VAL  :     L51  VAL    RMS: 0.446
    L52  SER  :     L52  SER    RMS: 0.766
    L53  ASN  :     L53  ASN    RMS: 0.755
    L54  ARG  :     L54  ARG    RMS: 1.860
    L55  PHE  :     L55  PHE    RMS: 1.145
    L56  SER  :     L56  SER    RMS: 1.006
    L57  GLY  :     L57  GLY    RMS: 0.882
    L58  VAL  :     L58  VAL    RMS: 0.738
    L59  PRO  :     L59  PRO    RMS: 1.342
    L60  ASP  :     L60  ASP    RMS: 1.664
    L61  ARG  :     L61  ARG    RMS: 1.517
    L62  PHE  :     L62  PHE    RMS: 1.968
    L63  SER  :     L63  SER    RMS: 1.140
    L64  GLY  :     L64  GLY    RMS: 0.645
    L65  SER  :     L65  SER    RMS: 0.702
    L66  GLY  :     L66  GLY    RMS: 0.592
    L67  SER  :     L67  SER    RMS: 1.796
    L68  GLY  :     L68  GLY    RMS: 0.762
    L69  THR  :     L69  THR    RMS: 0.577
    L70  ASP  :     L70  ASP    RMS: 1.191
    L71  PHE  :     L71  PHE    RMS: 0.419
    L72  THR  :     L72  THR    RMS: 0.624
    L73  LEU  :     L73  LEU    RMS: 0.551
    L74  LYS  :     L74  LYS    RMS: 1.160
    L75  ILE  :     L75  ILE    RMS: 0.830
    L76  SER  :     L76  SER    RMS: 1.558
    L77  ARG  :     L77  ARG    RMS: 1.680
    L78  VAL  :     L78  VAL    RMS: 1.167
    L79  GLU  :     L79  GLU    RMS: 1.213
    L80  ALA  :     L80  ALA    RMS: 1.405
    L81  GLU  :     L81  GLU    RMS: 2.841
    L82  ASP  :     L82  ASP    RMS: 0.850
    L83  LEU  :     L83  LEU    RMS: 1.942
    L84  GLY  :     L84  GLY    RMS: 0.221
    L85  VAL  :     L85  VAL    RMS: 0.235
    L86  TYR  :     L86  TYR    RMS: 0.424
    L87  PHE  :     L87  PHE    RMS: 0.146
    L88  CYS  :     L88  CYS    RMS: 0.223
    L89  SER  :     L89  SER    RMS: 0.455
    L90  GLN  :     L90  GLN    RMS: 1.001
    L91  SER  :     L91  SER    RMS: 0.955
    L92  THR  :     L92  THR    RMS: 0.872
    L93  HIS  :     L93  HIS    RMS: 1.752
    L94  VAL  :     L94  VAL    RMS: 3.278
    L95  PRO  :     L95  PRO    RMS: 1.977
    L96  PRO  :     L96  PRO    RMS: 0.617
    L97  THR  :     L97  THR    RMS: 0.427
    L98  PHE  :     L98  PHE    RMS: 1.511
    L99  GLY  :     L99  GLY    RMS: 0.313
   L100  GLY  :    L100  GLY    RMS: 0.380
   L101  GLY  :    L101  GLY    RMS: 0.345
   L102  THR  :    L102  THR    RMS: 0.280
   L103  LYS  :    L103  LYS    RMS: 0.697
   L104  LEU  :    L104  LEU    RMS: 0.463
   L105  GLU  :    L105  GLU    RMS: 1.541
   L106  ILE  :    L106  ILE    RMS: 2.797
   L107  LYS  :    L107  LYS    RMS: 2.179
   L108  ARG  :    L108  ARG    RMS: 7.182
   L109  THR  :    L109  THR    RMS: 7.879
   L110  VAL  :    L110  VAL    RMS: 10.810
     H1  GLN  :      H1  GLN    RMS: 2.210
     H2  VAL  :      H2  VAL    RMS: 0.835
     H3  GLN  :      H3  GLN    RMS: 1.816
     H4  LEU  :      H4  LEU    RMS: 0.905
     H5  LEU  :      H5  LEU    RMS: 1.661
     H6  GLU  :      H6  GLU    RMS: 1.305
     H7  SER  :      H7  SER    RMS: 0.671
     H8  GLY  :      H8  GLY    RMS: 0.883
     H9  PRO  :      H9  PRO    RMS: 0.844
    H10  GLU  :     H10  GLU    RMS: 0.890
    H11  LEU  :     H11  LEU    RMS: 1.119
    H12  LYS  :     H12  LYS    RMS: 1.085
    H13  LYS  :     H13  LYS    RMS: 1.207
    H14  PRO  :     H14  PRO    RMS: 0.834
    H15  GLY  :     H15  GLY    RMS: 0.734
    H16  GLU  :     H16  GLU    RMS: 0.739
    H17  THR  :     H17  THR    RMS: 0.737
    H18  VAL  :     H18  VAL    RMS: 0.530
    H19  LYS  :     H19  LYS    RMS: 0.846
    H20  ILE  :     H20  ILE    RMS: 0.654
    H21  SER  :     H21  SER    RMS: 0.617
    H22  CYS  :     H22  CYS    RMS: 0.361
    H23  LYS  :     H23  LYS    RMS: 1.176
    H24  ALA  :     H24  ALA    RMS: 0.391
    H25  SER  :     H25  SER    RMS: 1.105
    H26  GLY  :     H26  GLY    RMS: 2.198
    H27  TYR  :     H27  TYR    RMS: 2.381
    H28  THR  :     H28  THR    RMS: 2.999
    H29  PHE  :     H29  PHE    RMS: 1.176
    H30  THR  :     H30  THR    RMS: 1.689
    H31  ASN  :     H31  ASN    RMS: 1.855
    H32  TYR  :     H32  TYR    RMS: 1.076
    H33  GLY  :     H33  GLY    RMS: 0.804
    H34  MET  :     H34  MET    RMS: 1.419
    H35  ASN  :     H35  ASN    RMS: 0.535
    H36  TRP  :     H36  TRP    RMS: 0.510
    H37  VAL  :     H37  VAL    RMS: 0.546
    H38  LYS  :     H38  LYS    RMS: 0.720
    H39  GLN  :     H39  GLN    RMS: 0.631
    H40  ALA  :     H40  ALA    RMS: 0.902
    H41  PRO  :     H41  PRO    RMS: 1.424
    H42  GLY  :     H42  GLY    RMS: 1.309
    H43  LYS  :     H43  LYS    RMS: 1.232
    H44  GLY  :     H44  GLY    RMS: 1.119
    H45  LEU  :     H45  LEU    RMS: 0.461
    H46  LYS  :     H46  LYS    RMS: 0.941
    H47  TRP  :     H47  TRP    RMS: 0.788
    H48  MET  :     H48  MET    RMS: 0.373
    H49  GLY  :     H49  GLY    RMS: 0.296
    H50  TRP  :     H50  TRP    RMS: 0.642
    H51  ILE  :     H51  ILE    RMS: 0.746
    H52  ASN  :     H52  ASN    RMS: 0.907
    H52A THR  :     H52A THR    RMS: 1.159
    H53  TYR  :     H53  TYR    RMS: 2.138
    H54  THR  :     H54  THR    RMS: 1.578
    H55  GLY  :     H55  GLY    RMS: 0.935
    H56  GLU  :     H56  GLU    RMS: 1.632
    H57  PRO  :     H57  PRO    RMS: 0.754
    H58  THR  :     H58  THR    RMS: 0.915
    H59  TYR  :     H59  TYR    RMS: 0.782
    H60  ALA  :     H60  ALA    RMS: 0.648
    H61  ASP  :     H61  ASP    RMS: 1.541
    H62  ASP  :     H62  ASP    RMS: 1.458
    H63  PHE  :     H63  PHE    RMS: 1.626
    H64  LYS  :     H64  LYS    RMS: 1.247
    H65  GLY  :     H65  GLY    RMS: 0.853
    H66  ARG  :     H66  ARG    RMS: 0.839
    H67  PHE  :     H67  PHE    RMS: 0.710
    H68  ALA  :     H68  ALA    RMS: 0.596
    H69  PHE  :     H69  PHE    RMS: 0.435
    H70  SER  :     H70  SER    RMS: 0.671
    H71  LEU  :     H71  LEU    RMS: 0.693
    H72  GLU  :     H72  GLU    RMS: 1.283
    H73  THR  :     H73  THR    RMS: 1.577
    H74  SER  :     H74  SER    RMS: 0.946
    H75  ALA  :     H75  ALA    RMS: 0.520
    H76  SER  :     H76  SER    RMS: 0.866
    H77  THR  :     H77  THR    RMS: 0.587
    H78  ALA  :     H78  ALA    RMS: 0.490
    H79  TYR  :     H79  TYR    RMS: 1.538
    H80  LEU  :     H80  LEU    RMS: 0.467
    H81  GLN  :     H81  GLN    RMS: 0.986
    H82  ILE  :     H82  ILE    RMS: 0.630
    H82A ASN  :     H82A ASN    RMS: 1.406
    H82B ASN  :     H82B ASN    RMS: 2.560
    H82C LEU  :     H82C LEU    RMS: 1.070
    H83  LYS  :     H83  LYS    RMS: 1.533
    H84  ASN  :     H84  ASN    RMS: 0.895
    H85  GLU  :     H85  GLU    RMS: 1.166
    H86  ASP  :     H86  ASP    RMS: 0.757
    H87  THR  :     H87  THR    RMS: 0.783
    H88  ALA  :     H88  ALA    RMS: 0.678
    H89  THR  :     H89  THR    RMS: 0.659
    H90  TYR  :     H90  TYR    RMS: 0.757
    H91  PHE  :     H91  PHE    RMS: 0.376
    H92  CYS  :     H92  CYS    RMS: 0.434
    H93  VAL  :     H93  VAL    RMS: 0.395
    H94  GLN  :     H94  GLN    RMS: 0.670
    H95  ALA  :     H95  ALA    RMS: 0.532
    H96  GLU  :     H96  GLU    RMS: 0.677
    H97  ARG  :     H97  ARG    RMS: 1.399
    H98  LEU  :     H98  LEU    RMS: 0.363
    H99  ARG  :     H99  ARG    RMS: 0.434
   H100  ARG  :    H100  ARG    RMS: 1.159
   H100A THR  :    H100A THR    RMS: 0.567
   H100B PHE  :    H100B PHE    RMS: 0.544
   H101  ASP  :    H101  ASP    RMS: 1.109
   H102  TYR  :    H102  TYR    RMS: 0.689
   H103  TRP  :    H103  TRP    RMS: 0.700
   H104  GLY  :    H104  GLY    RMS: 0.584
   H105  ALA  :    H105  ALA    RMS: 0.612
   H106  GLY  :    H106  GLY    RMS: 0.416
   H107  THR  :    H107  THR    RMS: 0.476
   H108  THR  :    H108  THR    RMS: 0.938
   H109  VAL  :    H109  VAL    RMS: 0.706
   H110  THR  :    H110  THR    RMS: 0.851
   H111  VAL  :    H111  VAL    RMS: 0.765
   H112  SER  :    H112  SER    RMS: 0.756
   H113  SER  :    H113  SER    RMS: 1.771
   RMS: 1.551
   #global Ca-Ca
   Warning==> Structures have not yet been fitted.
   Fitting structures...
   RMS: 1.216
   #global Ca-Ca by residue
     L1  GLU  :      L1  GLU    RMS: 2.978
     L2  LEU  :      L2  LEU    RMS: 1.232
     L3  VAL  :      L3  VAL    RMS: 0.789
     L4  MET  :      L4  MET    RMS: 0.549
     L5  THR  :      L5  THR    RMS: 0.338
     L6  GLN  :      L6  GLN    RMS: 0.210
     L7  THR  :      L7  THR    RMS: 1.095
     L8  PRO  :      L8  PRO    RMS: 0.975
     L9  LEU  :      L9  LEU    RMS: 0.544
    L10  SER  :     L10  SER    RMS: 0.395
    L11  LEU  :     L11  LEU    RMS: 0.215
    L12  PRO  :     L12  PRO    RMS: 0.329
    L13  VAL  :     L13  VAL    RMS: 0.953
    L14  SER  :     L14  SER    RMS: 1.035
    L15  LEU  :     L15  LEU    RMS: 1.541
    L16  GLY  :     L16  GLY    RMS: 1.814
    L17  ASP  :     L17  ASP    RMS: 1.427
    L18  GLN  :     L18  GLN    RMS: 1.145
    L19  ALA  :     L19  ALA    RMS: 0.755
    L20  SER  :     L20  SER    RMS: 0.767
    L21  ILE  :     L21  ILE    RMS: 0.628
    L22  SER  :     L22  SER    RMS: 0.653
    L23  CYS  :     L23  CYS    RMS: 0.437
    L24  ARG  :     L24  ARG    RMS: 0.449
    L25  SER  :     L25  SER    RMS: 0.451
    L26  SER  :     L26  SER    RMS: 0.524
    L27  GLN  :     L27  GLN    RMS: 0.468
    L28  SER  :     L28  SER    RMS: 0.384
    L29  LEU  :     L29  LEU    RMS: 0.148
    L30  VAL  :     L30  VAL    RMS: 0.181
    L30A HIS  :     L30A HIS    RMS: 0.339
    L30B SER  :     L30B SER    RMS: 0.621
    L30C ASN  :     L30C ASN    RMS: 0.708
    L30D GLY  :     L30D GLY    RMS: 0.629
    L30E ASN  :     L30E ASN    RMS: 0.662
    L31  THR  :     L31  THR    RMS: 0.160
    L32  TYR  :     L32  TYR    RMS: 0.363
    L33  LEU  :     L33  LEU    RMS: 0.248
    L34  HIS  :     L34  HIS    RMS: 0.210
    L35  TRP  :     L35  TRP    RMS: 0.284
    L36  TYR  :     L36  TYR    RMS: 0.318
    L37  LEU  :     L37  LEU    RMS: 0.284
    L38  GLN  :     L38  GLN    RMS: 0.216
    L39  LYS  :     L39  LYS    RMS: 0.899
    L40  PRO  :     L40  PRO    RMS: 1.412
    L41  GLY  :     L41  GLY    RMS: 1.198
    L42  GLN  :     L42  GLN    RMS: 0.612
    L43  SER  :     L43  SER    RMS: 0.584
    L44  PRO  :     L44  PRO    RMS: 0.477
    L45  LYS  :     L45  LYS    RMS: 0.630
    L46  LEU  :     L46  LEU    RMS: 0.540
    L47  LEU  :     L47  LEU    RMS: 0.175
    L48  ILE  :     L48  ILE    RMS: 0.039
    L49  TYR  :     L49  TYR    RMS: 0.137
    L50  LYS  :     L50  LYS    RMS: 0.161
    L51  VAL  :     L51  VAL    RMS: 0.329
    L52  SER  :     L52  SER    RMS: 0.588
    L53  ASN  :     L53  ASN    RMS: 0.612
    L54  ARG  :     L54  ARG    RMS: 0.613
    L55  PHE  :     L55  PHE    RMS: 0.713
    L56  SER  :     L56  SER    RMS: 0.928
    L57  GLY  :     L57  GLY    RMS: 0.990
    L58  VAL  :     L58  VAL    RMS: 0.733
    L59  PRO  :     L59  PRO    RMS: 1.147
    L60  ASP  :     L60  ASP    RMS: 1.196
    L61  ARG  :     L61  ARG    RMS: 0.907
    L62  PHE  :     L62  PHE    RMS: 1.053
    L63  SER  :     L63  SER    RMS: 0.783
    L64  GLY  :     L64  GLY    RMS: 0.603
    L65  SER  :     L65  SER    RMS: 0.464
    L66  GLY  :     L66  GLY    RMS: 0.556
    L67  SER  :     L67  SER    RMS: 0.942
    L68  GLY  :     L68  GLY    RMS: 0.746
    L69  THR  :     L69  THR    RMS: 0.477
    L70  ASP  :     L70  ASP    RMS: 0.436
    L71  PHE  :     L71  PHE    RMS: 0.466
    L72  THR  :     L72  THR    RMS: 0.507
    L73  LEU  :     L73  LEU    RMS: 0.495
    L74  LYS  :     L74  LYS    RMS: 0.708
    L75  ILE  :     L75  ILE    RMS: 0.719
    L76  SER  :     L76  SER    RMS: 1.088
    L77  ARG  :     L77  ARG    RMS: 1.087
    L78  VAL  :     L78  VAL    RMS: 1.057
    L79  GLU  :     L79  GLU    RMS: 1.214
    L80  ALA  :     L80  ALA    RMS: 1.362
    L81  GLU  :     L81  GLU    RMS: 1.258
    L82  ASP  :     L82  ASP    RMS: 0.766
    L83  LEU  :     L83  LEU    RMS: 0.281
    L84  GLY  :     L84  GLY    RMS: 0.273
    L85  VAL  :     L85  VAL    RMS: 0.167
    L86  TYR  :     L86  TYR    RMS: 0.199
    L87  PHE  :     L87  PHE    RMS: 0.080
    L88  CYS  :     L88  CYS    RMS: 0.146
    L89  SER  :     L89  SER    RMS: 0.260
    L90  GLN  :     L90  GLN    RMS: 0.346
    L91  SER  :     L91  SER    RMS: 0.639
    L92  THR  :     L92  THR    RMS: 0.694
    L93  HIS  :     L93  HIS    RMS: 1.442
    L94  VAL  :     L94  VAL    RMS: 2.251
    L95  PRO  :     L95  PRO    RMS: 1.337
    L96  PRO  :     L96  PRO    RMS: 0.529
    L97  THR  :     L97  THR    RMS: 0.390
    L98  PHE  :     L98  PHE    RMS: 0.234
    L99  GLY  :     L99  GLY    RMS: 0.315
   L100  GLY  :    L100  GLY    RMS: 0.424
   L101  GLY  :    L101  GLY    RMS: 0.149
   L102  THR  :    L102  THR    RMS: 0.125
   L103  LYS  :    L103  LYS    RMS: 0.084
   L104  LEU  :    L104  LEU    RMS: 0.311
   L105  GLU  :    L105  GLU    RMS: 0.583
   L106  ILE  :    L106  ILE    RMS: 1.925
   L107  LYS  :    L107  LYS    RMS: 2.434
   L108  ARG  :    L108  ARG    RMS: 5.159
   L109  THR  :    L109  THR    RMS: 7.485
   L110  VAL  :    L110  VAL    RMS: 10.023
     H1  GLN  :      H1  GLN    RMS: 1.415
     H2  VAL  :      H2  VAL    RMS: 0.571
     H3  GLN  :      H3  GLN    RMS: 0.622
     H4  LEU  :      H4  LEU    RMS: 0.773
     H5  LEU  :      H5  LEU    RMS: 1.116
     H6  GLU  :      H6  GLU    RMS: 0.098
     H7  SER  :      H7  SER    RMS: 0.632
     H8  GLY  :      H8  GLY    RMS: 0.825
     H9  PRO  :      H9  PRO    RMS: 0.877
    H10  GLU  :     H10  GLU    RMS: 0.831
    H11  LEU  :     H11  LEU    RMS: 1.056
    H12  LYS  :     H12  LYS    RMS: 1.061
    H13  LYS  :     H13  LYS    RMS: 0.935
    H14  PRO  :     H14  PRO    RMS: 0.819
    H15  GLY  :     H15  GLY    RMS: 0.791
    H16  GLU  :     H16  GLU    RMS: 0.697
    H17  THR  :     H17  THR    RMS: 0.675
    H18  VAL  :     H18  VAL    RMS: 0.556
    H19  LYS  :     H19  LYS    RMS: 0.678
    H20  ILE  :     H20  ILE    RMS: 0.518
    H21  SER  :     H21  SER    RMS: 0.485
    H22  CYS  :     H22  CYS    RMS: 0.326
    H23  LYS  :     H23  LYS    RMS: 0.353
    H24  ALA  :     H24  ALA    RMS: 0.304
    H25  SER  :     H25  SER    RMS: 1.003
    H26  GLY  :     H26  GLY    RMS: 2.155
    H27  TYR  :     H27  TYR    RMS: 1.541
    H28  THR  :     H28  THR    RMS: 2.799
    H29  PHE  :     H29  PHE    RMS: 1.294
    H30  THR  :     H30  THR    RMS: 1.470
    H31  ASN  :     H31  ASN    RMS: 1.768
    H32  TYR  :     H32  TYR    RMS: 1.037
    H33  GLY  :     H33  GLY    RMS: 0.961
    H34  MET  :     H34  MET    RMS: 0.765
    H35  ASN  :     H35  ASN    RMS: 0.342
    H36  TRP  :     H36  TRP    RMS: 0.262
    H37  VAL  :     H37  VAL    RMS: 0.206
    H38  LYS  :     H38  LYS    RMS: 0.513
    H39  GLN  :     H39  GLN    RMS: 0.457
    H40  ALA  :     H40  ALA    RMS: 0.859
    H41  PRO  :     H41  PRO    RMS: 1.265
    H42  GLY  :     H42  GLY    RMS: 1.369
    H43  LYS  :     H43  LYS    RMS: 0.791
    H44  GLY  :     H44  GLY    RMS: 0.474
    H45  LEU  :     H45  LEU    RMS: 0.395
    H46  LYS  :     H46  LYS    RMS: 0.332
    H47  TRP  :     H47  TRP    RMS: 0.159
    H48  MET  :     H48  MET    RMS: 0.378
    H49  GLY  :     H49  GLY    RMS: 0.366
    H50  TRP  :     H50  TRP    RMS: 0.354
    H51  ILE  :     H51  ILE    RMS: 0.779
    H52  ASN  :     H52  ASN    RMS: 0.963
    H52A THR  :     H52A THR    RMS: 1.144
    H53  TYR  :     H53  TYR    RMS: 1.522
    H54  THR  :     H54  THR    RMS: 1.168
    H55  GLY  :     H55  GLY    RMS: 0.862
    H56  GLU  :     H56  GLU    RMS: 0.984
    H57  PRO  :     H57  PRO    RMS: 0.790
    H58  THR  :     H58  THR    RMS: 0.886
    H59  TYR  :     H59  TYR    RMS: 0.715
    H60  ALA  :     H60  ALA    RMS: 0.619
    H61  ASP  :     H61  ASP    RMS: 1.104
    H62  ASP  :     H62  ASP    RMS: 1.010
    H63  PHE  :     H63  PHE    RMS: 0.736
    H64  LYS  :     H64  LYS    RMS: 0.693
    H65  GLY  :     H65  GLY    RMS: 0.871
    H66  ARG  :     H66  ARG    RMS: 0.720
    H67  PHE  :     H67  PHE    RMS: 0.479
    H68  ALA  :     H68  ALA    RMS: 0.584
    H69  PHE  :     H69  PHE    RMS: 0.447
    H70  SER  :     H70  SER    RMS: 0.570
    H71  LEU  :     H71  LEU    RMS: 0.464
    H72  GLU  :     H72  GLU    RMS: 0.582
    H73  THR  :     H73  THR    RMS: 1.127
    H74  SER  :     H74  SER    RMS: 0.877
    H75  ALA  :     H75  ALA    RMS: 0.423
    H76  SER  :     H76  SER    RMS: 0.380
    H77  THR  :     H77  THR    RMS: 0.432
    H78  ALA  :     H78  ALA    RMS: 0.469
    H79  TYR  :     H79  TYR    RMS: 0.574
    H80  LEU  :     H80  LEU    RMS: 0.484
    H81  GLN  :     H81  GLN    RMS: 0.657
    H82  ILE  :     H82  ILE    RMS: 0.635
    H82A ASN  :     H82A ASN    RMS: 0.778
    H82B ASN  :     H82B ASN    RMS: 1.082
    H82C LEU  :     H82C LEU    RMS: 0.979
    H83  LYS  :     H83  LYS    RMS: 0.729
    H84  ASN  :     H84  ASN    RMS: 0.874
    H85  GLU  :     H85  GLU    RMS: 1.209
    H86  ASP  :     H86  ASP    RMS: 0.722
    H87  THR  :     H87  THR    RMS: 0.736
    H88  ALA  :     H88  ALA    RMS: 0.613
    H89  THR  :     H89  THR    RMS: 0.601
    H90  TYR  :     H90  TYR    RMS: 0.561
    H91  PHE  :     H91  PHE    RMS: 0.373
    H92  CYS  :     H92  CYS    RMS: 0.323
    H93  VAL  :     H93  VAL    RMS: 0.356
    H94  GLN  :     H94  GLN    RMS: 0.484
    H95  ALA  :     H95  ALA    RMS: 0.552
    H96  GLU  :     H96  GLU    RMS: 0.640
    H97  ARG  :     H97  ARG    RMS: 0.617
    H98  LEU  :     H98  LEU    RMS: 0.123
    H99  ARG  :     H99  ARG    RMS: 0.216
   H100  ARG  :    H100  ARG    RMS: 0.389
   H100A THR  :    H100A THR    RMS: 0.608
   H100B PHE  :    H100B PHE    RMS: 0.426
   H101  ASP  :    H101  ASP    RMS: 0.437
   H102  TYR  :    H102  TYR    RMS: 0.628
   H103  TRP  :    H103  TRP    RMS: 0.612
   H104  GLY  :    H104  GLY    RMS: 0.605
   H105  ALA  :    H105  ALA    RMS: 0.589
   H106  GLY  :    H106  GLY    RMS: 0.361
   H107  THR  :    H107  THR    RMS: 0.468
   H108  THR  :    H108  THR    RMS: 0.700
   H109  VAL  :    H109  VAL    RMS: 0.696
   H110  THR  :    H110  THR    RMS: 0.759
   H111  VAL  :    H111  VAL    RMS: 0.806
   H112  SER  :    H112  SER    RMS: 0.504
   H113  SER  :    H113  SER    RMS: 0.809
   RMS: 1.216
'''
