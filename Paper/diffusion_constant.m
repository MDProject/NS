% Given mean square distance 
% 0-0-0.5   x
Eins = [0.000000	0.000000
0.020000	0.000688
0.040000	0.002609
0.060000	0.005431
0.080000	0.008817
0.100000	0.012507
0.120000	0.016341
0.140000	0.020235
0.160000	0.024155
0.180000	0.028095
0.200000	0.032058
0.220000	0.036048
0.240000	0.040067
0.260000	0.044112
0.280000	0.048179
0.300000	0.052263
0.320000	0.056360
0.340000	0.060467
0.360000	0.064584
0.380000	0.068712
0.400000	0.072856
0.420000	0.077018
0.440000	0.081198
0.460000	0.085396
0.480000	0.089608
0.500000	0.093833
0.520000	0.098068
0.540000	0.102311
0.560000	0.106562
0.580000	0.110820
0.600000	0.115082
0.620000	0.119348
0.640000	0.123615
0.660000	0.127883
0.680000	0.132152
0.700000	0.136424
0.720000	0.140698
0.740000	0.144978
0.760000	0.149263
0.780000	0.153554
0.800000	0.157849
0.820000	0.162146
0.840000	0.166444
0.860000	0.170741
0.880000	0.175036
0.900000	0.179334
0.920000	0.183636
0.940000	0.187946
0.960000	0.192267
0.980000	0.196604
1.000000	0.200959
1.020000	0.205335
1.040000	0.209733
1.060000	0.214151
1.080000	0.218586
1.100000	0.223028
1.120000	0.227471
1.140000	0.231908
1.160000	0.236334
1.180000	0.240748
1.200000	0.245147
1.220000	0.249530
1.240000	0.253898
1.260000	0.258251
1.280000	0.262594
1.300000	0.266932
1.320000	0.271269
1.340000	0.275611
1.360000	0.279960
1.380000	0.284320
1.400000	0.288685
1.420000	0.293048
1.440000	0.297400
1.460000	0.301736
1.480000	0.306053
1.500000	0.310351
1.520000	0.314633
1.540000	0.318898
1.560000	0.323150
1.580000	0.327387];
time = Eins(:,1);
Diff = Eins(:,2);
t = time(floor(length(time)/10*7):length(time));
D = Diff(floor(length(time)/10*7):length(time)); % 0.218
figure('visible','on');
plot(time,Diff,'or','linewidth',2.7);
hold on;
plot(time,0.218*(time-t(1))+D(1),'-b','linewidth',2.2);
set(gca,'fontsize', 60);xlabel({'$\Delta t\ (\sqrt{m\sigma^2/\varepsilon})$'},'fontsize',50,'Interpreter','latex');ylabel({'$1/2<[x(t)-x(0)]^2> (\sigma^2)$'},'fontsize',50,'Interpreter','latex');

% 0-2-0.5  x
Eins = [0.000000	0.000000
0.020000	0.000689
0.040000	0.002614
0.060000	0.005458
0.080000	0.008890
0.100000	0.012659
0.120000	0.016605
0.140000	0.020643
0.160000	0.024731
0.180000	0.028851
0.200000	0.032995
0.220000	0.037158
0.240000	0.041337
0.260000	0.045529
0.280000	0.049737
0.300000	0.053964
0.320000	0.058216
0.340000	0.062501
0.360000	0.066823
0.380000	0.071188
0.400000	0.075597
0.420000	0.080048
0.440000	0.084538
0.460000	0.089063
0.480000	0.093618
0.500000	0.098205
0.520000	0.102826
0.540000	0.107482
0.560000	0.112172
0.580000	0.116889
0.600000	0.121626
0.620000	0.126375
0.640000	0.131131
0.660000	0.135891
0.680000	0.140652
0.700000	0.145413
0.720000	0.150178
0.740000	0.154945
0.760000	0.159709
0.780000	0.164460
0.800000	0.169189
0.820000	0.173892
0.840000	0.178563
0.860000	0.183202
0.880000	0.187805
0.900000	0.192373
0.920000	0.196907
0.940000	0.201410
0.960000	0.205885
0.980000	0.210331
1.000000	0.214749
1.020000	0.219140
1.040000	0.223504
1.060000	0.227836
1.080000	0.232135
1.100000	0.236398
1.120000	0.240631
1.140000	0.244840
1.160000	0.249034
1.180000	0.253220
1.200000	0.257401
1.220000	0.261579
1.240000	0.265754
1.260000	0.269923
1.280000	0.274089
1.300000	0.278253
1.320000	0.282421
1.340000	0.286595
1.360000	0.290777
1.380000	0.294966
1.400000	0.299156
1.420000	0.303336
1.440000	0.307495
1.460000	0.311628
1.480000	0.315733
1.500000	0.319811
1.520000	0.323871
1.540000	0.327925
1.560000	0.331981
1.580000	0.336048];
time = Eins(:,1);
Diff = Eins(:,2);
t = time(floor(length(time)/10*7):length(time));
D = Diff(floor(length(time)/10*7):length(time)); % 0.208


% 0-4-0.5  x
Eins = [0.000000	0.000000
0.020000	0.000700
0.040000	0.002650
0.060000	0.005508
0.080000	0.008929
0.100000	0.012652
0.120000	0.016516
0.140000	0.020438
0.160000	0.024384
0.180000	0.028347
0.200000	0.032335
0.220000	0.036355
0.240000	0.040416
0.260000	0.044522
0.280000	0.048671
0.300000	0.052856
0.320000	0.057072
0.340000	0.061313
0.360000	0.065578
0.380000	0.069875
0.400000	0.074211
0.420000	0.078591
0.440000	0.083014
0.460000	0.087476
0.480000	0.091973
0.500000	0.096498
0.520000	0.101050
0.540000	0.105623
0.560000	0.110211
0.580000	0.114808
0.600000	0.119409
0.620000	0.124013
0.640000	0.128620
0.660000	0.133234
0.680000	0.137860
0.700000	0.142501
0.720000	0.147162
0.740000	0.151845
0.760000	0.156550
0.780000	0.161276
0.800000	0.166018
0.820000	0.170773
0.840000	0.175536
0.860000	0.180301
0.880000	0.185067
0.900000	0.189834
0.920000	0.194604
0.940000	0.199376
0.960000	0.204151
0.980000	0.208925
1.000000	0.213695
1.020000	0.218458
1.040000	0.223211
1.060000	0.227957
1.080000	0.232699
1.100000	0.237442
1.120000	0.242194
1.140000	0.246959
1.160000	0.251738
1.180000	0.256528
1.200000	0.261324
1.220000	0.266122
1.240000	0.270918
1.260000	0.275712
1.280000	0.280501
1.300000	0.285283
1.320000	0.290056
1.340000	0.294819
1.360000	0.299574
1.380000	0.304322
1.400000	0.309067
1.420000	0.313809
1.440000	0.318549
1.460000	0.323289
1.480000	0.328032
1.500000	0.332776
1.520000	0.337524
1.540000	0.342273
1.560000	0.347025
1.580000	0.351781];
time = Eins(:,1);
Diff = Eins(:,2);
t = time(floor(length(time)/10*7):length(time));
D = Diff(floor(length(time)/10*7):length(time)); % 0.238


% 0-6-0.5  x
Eins = [0.000000	0.000000
0.020000	0.000699
0.040000	0.002656
0.060000	0.005546
0.080000	0.009035
0.100000	0.012864
0.120000	0.016863
0.140000	0.020942
0.160000	0.025058
0.180000	0.029199
0.200000	0.033363
0.220000	0.037556
0.240000	0.041786
0.260000	0.046057
0.280000	0.050377
0.300000	0.054748
0.320000	0.059173
0.340000	0.063654
0.360000	0.068189
0.380000	0.072778
0.400000	0.077416
0.420000	0.082100
0.440000	0.086821
0.460000	0.091574
0.480000	0.096356
0.500000	0.101167
0.520000	0.106010
0.540000	0.110890
0.560000	0.115811
0.580000	0.120772
0.600000	0.125775
0.620000	0.130819
0.640000	0.135904
0.660000	0.141023
0.680000	0.146170
0.700000	0.151333
0.720000	0.156503
0.740000	0.161674
0.760000	0.166843
0.780000	0.172013
0.800000	0.177190
0.820000	0.182377
0.840000	0.187572
0.860000	0.192769
0.880000	0.197956
0.900000	0.203123
0.920000	0.208263
0.940000	0.213375
0.960000	0.218459
0.980000	0.223521
1.000000	0.228566
1.020000	0.233598
1.040000	0.238614
1.060000	0.243613
1.080000	0.248591
1.100000	0.253549
1.120000	0.258489
1.140000	0.263416
1.160000	0.268333
1.180000	0.273240
1.200000	0.278132
1.220000	0.283005
1.240000	0.287853
1.260000	0.292674
1.280000	0.297470
1.300000	0.302242
1.320000	0.306992
1.340000	0.311722
1.360000	0.316433
1.380000	0.321129
1.400000	0.325814
1.420000	0.330495
1.440000	0.335179
1.460000	0.339870
1.480000	0.344571
1.500000	0.349284
1.520000	0.354011
1.540000	0.358757
1.560000	0.363526
1.580000	0.368320];
time = Eins(:,1);
Diff = Eins(:,2);
t = time(floor(length(time)/10*7):length(time));
D = Diff(floor(length(time)/10*7):length(time)); % 0.237


% 0-8-0.5  x
Eins = [0.000000	0.000000
0.020000	0.000694
0.040000	0.002627
0.060000	0.005466
0.080000	0.008864
0.100000	0.012553
0.120000	0.016367
0.140000	0.020222
0.160000	0.024083
0.180000	0.027942
0.200000	0.031802
0.220000	0.035670
0.240000	0.039552
0.260000	0.043454
0.280000	0.047385
0.300000	0.051354
0.320000	0.055369
0.340000	0.059435
0.360000	0.063552
0.380000	0.067716
0.400000	0.071918
0.420000	0.076151
0.440000	0.080406
0.460000	0.084682
0.480000	0.088976
0.500000	0.093286
0.520000	0.097607
0.540000	0.101932
0.560000	0.106256
0.580000	0.110576
0.600000	0.114896
0.620000	0.119220
0.640000	0.123557
0.660000	0.127911
0.680000	0.132286
0.700000	0.136684
0.720000	0.141100
0.740000	0.145531
0.760000	0.149973
0.780000	0.154427
0.800000	0.158896
0.820000	0.163381
0.840000	0.167883
0.860000	0.172399
0.880000	0.176928
0.900000	0.181464
0.920000	0.186005
0.940000	0.190550
0.960000	0.195101
0.980000	0.199660
1.000000	0.204232
1.020000	0.208818
1.040000	0.213417
1.060000	0.218024
1.080000	0.222631
1.100000	0.227231
1.120000	0.231815
1.140000	0.236377
1.160000	0.240912
1.180000	0.245421
1.200000	0.249910
1.220000	0.254386
1.240000	0.258854
1.260000	0.263315
1.280000	0.267768
1.300000	0.272213
1.320000	0.276650
1.340000	0.281083
1.360000	0.285515
1.380000	0.289947
1.400000	0.294383
1.420000	0.298823
1.440000	0.303273
1.460000	0.307732
1.480000	0.312200
1.500000	0.316672
1.520000	0.321141
1.540000	0.325603
1.560000	0.330061
1.580000	0.334518];
time = Eins(:,1);
Diff = Eins(:,2);
t = time(floor(length(time)/10*7):length(time));
D = Diff(floor(length(time)/10*7):length(time)); % 0.222



% 0-10-0.5  x
Eins = [0.000000	0.000000
0.020000	0.000683
0.040000	0.002587
0.060000	0.005381
0.080000	0.008726
0.100000	0.012365
0.120000	0.016135
0.140000	0.019955
0.160000	0.023791
0.180000	0.027638
0.200000	0.031503
0.220000	0.035396
0.240000	0.039322
0.260000	0.043283
0.280000	0.047281
0.300000	0.051313
0.320000	0.055380
0.340000	0.059480
0.360000	0.063612
0.380000	0.067771
0.400000	0.071953
0.420000	0.076153
0.440000	0.080368
0.460000	0.084597
0.480000	0.088838
0.500000	0.093095
0.520000	0.097368
0.540000	0.101659
0.560000	0.105965
0.580000	0.110282
0.600000	0.114610
0.620000	0.118947
0.640000	0.123298
0.660000	0.127664
0.680000	0.132050
0.700000	0.136457
0.720000	0.140885
0.740000	0.145333
0.760000	0.149794
0.780000	0.154262
0.800000	0.158729
0.820000	0.163188
0.840000	0.167636
0.860000	0.172075
0.880000	0.176509
0.900000	0.180946
0.920000	0.185389
0.940000	0.189839
0.960000	0.194297
0.980000	0.198761
1.000000	0.203231
1.020000	0.207707
1.040000	0.212191
1.060000	0.216681
1.080000	0.221178
1.100000	0.225684
1.120000	0.230199
1.140000	0.234726
1.160000	0.239265
1.180000	0.243817
1.200000	0.248378
1.220000	0.252941
1.240000	0.257501
1.260000	0.262054
1.280000	0.266596
1.300000	0.271124
1.320000	0.275637
1.340000	0.280135
1.360000	0.284619
1.380000	0.289093
1.400000	0.293559
1.420000	0.298019
1.440000	0.302473
1.460000	0.306923
1.480000	0.311367
1.500000	0.315805
1.520000	0.320231
1.540000	0.324645
1.560000	0.329042
1.580000	0.333422];
time = Eins(:,1);
Diff = Eins(:,2);
t = time(floor(length(time)/10*7):length(time));
D = Diff(floor(length(time)/10*7):length(time)); % 0.224


% 0-12-0.5  x
Eins = [0.000000	0.000000
0.020000	0.000670
0.040000	0.002540
0.060000	0.005289
0.080000	0.008586
0.100000	0.012175
0.120000	0.015896
0.140000	0.019666
0.160000	0.023453
0.180000	0.027251
0.200000	0.031064
0.220000	0.034901
0.240000	0.038764
0.260000	0.042657
0.280000	0.046584
0.300000	0.050548
0.320000	0.054548
0.340000	0.058582
0.360000	0.062643
0.380000	0.066725
0.400000	0.070821
0.420000	0.074927
0.440000	0.079039
0.460000	0.083156
0.480000	0.087280
0.500000	0.091409
0.520000	0.095547
0.540000	0.099697
0.560000	0.103864
0.580000	0.108051
0.600000	0.112262
0.620000	0.116497
0.640000	0.120757
0.660000	0.125037
0.680000	0.129336
0.700000	0.133649
0.720000	0.137975
0.740000	0.142307
0.760000	0.146638
0.780000	0.150962
0.800000	0.155278
0.820000	0.159585
0.840000	0.163888
0.860000	0.168192
0.880000	0.172504
0.900000	0.176830
0.920000	0.181176
0.940000	0.185543
0.960000	0.189934
0.980000	0.194346
1.000000	0.198774
1.020000	0.203211
1.040000	0.207646
1.060000	0.212067
1.080000	0.216468
1.100000	0.220845
1.120000	0.225202
1.140000	0.229539
1.160000	0.233858
1.180000	0.238161
1.200000	0.242450
1.220000	0.246731
1.240000	0.251010
1.260000	0.255293
1.280000	0.259585
1.300000	0.263886
1.320000	0.268194
1.340000	0.272509
1.360000	0.276830
1.380000	0.281161
1.400000	0.285503
1.420000	0.289858
1.440000	0.294231
1.460000	0.298621
1.480000	0.303031
1.500000	0.307461
1.520000	0.311908
1.540000	0.316368
1.560000	0.320837
1.580000	0.325310];
time = Eins(:,1);
Diff = Eins(:,2);
t = time(floor(length(time)/10*7):length(time));
D = Diff(floor(length(time)/10*7):length(time)); % 0.217


% 0-14-0.5  x
Eins = [0.000000	0.000000
0.020000	0.000687
0.040000	0.002604
0.060000	0.005422
0.080000	0.008805
0.100000	0.012496
0.120000	0.016333
0.140000	0.020229
0.160000	0.024149
0.180000	0.028084
0.200000	0.032035
0.220000	0.036009
0.240000	0.040013
0.260000	0.044055
0.280000	0.048145
0.300000	0.052285
0.320000	0.056472
0.340000	0.060699
0.360000	0.064960
0.380000	0.069249
0.400000	0.073565
0.420000	0.077910
0.440000	0.082287
0.460000	0.086698
0.480000	0.091143
0.500000	0.095621
0.520000	0.100127
0.540000	0.104659
0.560000	0.109215
0.580000	0.113792
0.600000	0.118387
0.620000	0.122997
0.640000	0.127620
0.660000	0.132252
0.680000	0.136892
0.700000	0.141537
0.720000	0.146185
0.740000	0.150835
0.760000	0.155484
0.780000	0.160131
0.800000	0.164773
0.820000	0.169413
0.840000	0.174051
0.860000	0.178690
0.880000	0.183331
0.900000	0.187977
0.920000	0.192635
0.940000	0.197307
0.960000	0.201997
0.980000	0.206705
1.000000	0.211428
1.020000	0.216162
1.040000	0.220904
1.060000	0.225650
1.080000	0.230401
1.100000	0.235158
1.120000	0.239926
1.140000	0.244707
1.160000	0.249502
1.180000	0.254311
1.200000	0.259131
1.220000	0.263960
1.240000	0.268787
1.260000	0.273605
1.280000	0.278406
1.300000	0.283183
1.320000	0.287938
1.340000	0.292671
1.360000	0.297385
1.380000	0.302082
1.400000	0.306762
1.420000	0.311427
1.440000	0.316079
1.460000	0.320723
1.480000	0.325364
1.500000	0.330007
1.520000	0.334655
1.540000	0.339307
1.560000	0.343960
1.580000	0.348614];
time = Eins(:,1);
Diff = Eins(:,2);
t = time(floor(length(time)/10*7):length(time));
D = Diff(floor(length(time)/10*7):length(time)); % 0.235




% 0-16-0.5  x
Eins = [0.000000	0.000000
0.020000	0.000690
0.040000	0.002618
0.060000	0.005454
0.080000	0.008860
0.100000	0.012571
0.120000	0.016422
0.140000	0.020323
0.160000	0.024237
0.180000	0.028155
0.200000	0.032077
0.220000	0.036007
0.240000	0.039951
0.260000	0.043909
0.280000	0.047883
0.300000	0.051867
0.320000	0.055859
0.340000	0.059855
0.360000	0.063857
0.380000	0.067866
0.400000	0.071886
0.420000	0.075916
0.440000	0.079955
0.460000	0.084001
0.480000	0.088054
0.500000	0.092113
0.520000	0.096176
0.540000	0.100236
0.560000	0.104290
0.580000	0.108333
0.600000	0.112363
0.620000	0.116380
0.640000	0.120387
0.660000	0.124388
0.680000	0.128386
0.700000	0.132382
0.720000	0.136374
0.740000	0.140360
0.760000	0.144340
0.780000	0.148316
0.800000	0.152294
0.820000	0.156280
0.840000	0.160277
0.860000	0.164291
0.880000	0.168321
0.900000	0.172373
0.920000	0.176451
0.940000	0.180558
0.960000	0.184697
0.980000	0.188866
1.000000	0.193060
1.020000	0.197269
1.040000	0.201489
1.060000	0.205713
1.080000	0.209941
1.100000	0.214171
1.120000	0.218404
1.140000	0.222638
1.160000	0.226877
1.180000	0.231121
1.200000	0.235373
1.220000	0.239632
1.240000	0.243898
1.260000	0.248171
1.280000	0.252446
1.300000	0.256720
1.320000	0.260990
1.340000	0.265254
1.360000	0.269511
1.380000	0.273760
1.400000	0.278003
1.420000	0.282244
1.440000	0.286484
1.460000	0.290724
1.480000	0.294963
1.500000	0.299198
1.520000	0.303430
1.540000	0.307656
1.560000	0.311875
1.580000	0.316084];
time = Eins(:,1);
Diff = Eins(:,2);
t = time(floor(length(time)/10*2):floor(length(time)/10*5));
D = Diff(floor(length(time)/10*2):floor(length(time)/10*5)); % 0.201
plot(time(1:floor(length(time)/10*5)),Diff(1:floor(length(time)/10*5)),'or','linewidth',2.2);
hold on;
plot(time(1:floor(length(time)/10*5)),0.201*(time(1:floor(length(time)/10*5))-t(1))+D(1),'-b','linewidth',2.2);
set(gca,'fontsize', 50);xlabel({'$Time\ interval\  \Delta t\ (in\ unit\ of\ \sqrt{\frac{m\sigma^2}{\varepsilon}})$'},'fontsize',50,'Interpreter','latex');ylabel({'$Half\ mean\ squared\ distance\ \frac{1}{2}<x(t)-x(0)>^2 (in\ unit\ of\ \sigma^2)$'},'fontsize',32,'Interpreter','latex');


% 0-17-0.5  x
clear all;
Eins = [0.000000	0.000000
0.020000	0.000706
0.040000	0.002673
0.060000	0.005558
0.080000	0.009007
0.100000	0.012754
0.120000	0.016628
0.140000	0.020542
0.160000	0.024458
0.180000	0.028366
0.200000	0.032270
0.220000	0.036180
0.240000	0.040101
0.260000	0.044035
0.280000	0.047975
0.300000	0.051912
0.320000	0.055837
0.340000	0.059747
0.360000	0.063644
0.380000	0.067532
0.400000	0.071419
0.420000	0.075307
0.440000	0.079195
0.460000	0.083086
0.480000	0.086982
0.500000	0.090884
0.520000	0.094793
0.540000	0.098707
0.560000	0.102625
0.580000	0.106542
0.600000	0.110448];
time = Eins(:,1);
Diff = Eins(:,2);
t = time(floor(length(time)/10*7):length(time));
D = Diff(floor(length(time)/10*7):length(time)); % 0.196
plot(time,Diff,'or','linewidth',2.2);
hold on;
plot(time,0.196*(time-t(1))+D(1),'-b','linewidth',2.2);
set(gca,'fontsize', 50);xlabel({'$Time\ interval\  \Delta t\ (in\ unit\ of\ \sqrt{\frac{m\sigma^2}{\varepsilon}})$'},'fontsize',50,'Interpreter','latex');ylabel({'$Half\ mean\ squared\ distance\ \frac{1}{2}<x(t)-x(0)>^2 (in\ unit\ of\ \sigma^2)$'},'fontsize',32,'Interpreter','latex');


% Kubo formula
% 0-0-0.4   x
clear all;
Kubo = [0.000000	0.495362
0.020000	0.439767
0.040000	0.317606
0.060000	0.194804
0.080000	0.103136
0.100000	0.047025
0.120000	0.018969
0.140000	0.008643
0.160000	0.007065
0.180000	0.008201
0.200000	0.009279
0.220000	0.009881
0.240000	0.010261
0.260000	0.010620
0.280000	0.011240
0.300000	0.012118
0.320000	0.012794
0.340000	0.012783
0.360000	0.011849
0.380000	0.010033
0.400000	0.007690
0.420000	0.005375
0.440000	0.003813
0.460000	0.003481
0.480000	0.004234
0.500000	0.005349
0.520000	0.005955
0.540000	0.005549
0.560000	0.004148
0.580000	0.002176
0.600000	0.000298
0.620000	-0.000864
0.640000	-0.001092
0.660000	-0.000716
0.680000	-0.000452
0.700000	-0.000579
0.720000	-0.000698
0.740000	-0.000465
0.760000	0.000128
0.780000	0.000788
0.800000	0.001457
0.820000	0.002270
0.840000	0.003161
0.860000	0.003881
0.880000	0.004158
0.900000	0.003726
0.920000	0.002825
0.940000	0.002104
0.960000	0.001862
0.980000	0.002043
1.000000	0.002453
1.020000	0.002757
1.040000	0.002660
1.060000	0.002271
1.080000	0.002106
1.100000	0.002480
1.120000	0.003250
1.140000	0.004015
1.160000	0.004457
1.180000	0.004444
1.200000	0.003942
1.220000	0.003098
1.240000	0.002287
1.260000	0.001709
1.280000	0.001334
1.300000	0.001109
1.320000	0.001035
1.340000	0.001063
1.360000	0.001104
1.380000	0.000955
1.400000	0.000760
1.420000	0.001032
1.440000	0.001856
1.460000	0.002634
1.480000	0.002787
1.500000	0.002333
1.520000	0.001562
1.540000	0.000614
1.560000	-0.000340
1.580000	-0.001165];
Diff = Kubo(:,2);
time = Kubo(:,1);
t = time(1:floor(length(time)/10*8));
D = Diff(1:floor(length(time)/10*8));
plot(t,D,'or','linewidth',3);
set(gca,'fontsize', 50);xlabel({'$\Delta t\ (\sqrt{m\sigma^2/\varepsilon})$'},'fontsize',40,'Interpreter','latex');ylabel({'$<\vec{v}(0)\cdot\vec{v}(t)>\ (\varepsilon/m$'},'fontsize',40,'Interpreter','latex');