function standard_R_RC_element_data()
  f = [0.03125, 0.0625, 0.125, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0, 256.0, 512.0, 1024.0, 2048.0, 4096.0, 8192.0, 16384.0, 32768.0,        
  65536.0, 131072.0, 262144.0]
  Z = [7.9999566250062575 - 0.014804038716500654im, 7.999826504539605 - 0.029607311013954685im, 7.999306090382245 - 0.05920849147321719im, 7.997225516509134 - 0.11836796402074179im, 7.988920507331971 - 0.2363445912343561im, 7.955974654576409 - 0.46958432924037374im, 7.828430016246774 - 0.9151288224796308im, 7.377764326740351 - 1.6603747779009874im, 6.187578798483693 - 2.423483564393849im, 4.526395789326643 - 2.3428639092377055im, 3.4920252494595356 - 1.5733429341387992im, 3.120896373931658 - 0.9653675926047738im, 2.983681817667524 - 0.7228105806992193im, 2.8210053520607663 - 0.785109509619991im, 2.4166123325165447 - 1.00398585317664im, 1.7535928232882099 - 1.0015305103631473im, 1.2624850542384358 - 0.6624837626405642im, 1.0727779748610742 - 0.3101732217968406im, 1.0187045081983228 - 0.04452524954644659im, 1.0047091282613076 + 0.208681260966955im, 1.0011793628768881 + 0.5675016337277559im, 1.0002949710575257 + 1.2102334554501843im, 1.0000737509150293 + 2.4581021305355573im, 1.000018438238243 + 4.935024388053915im]
  return f, Z
end

function extrapol_all_data()
  Z = Any[59.9998025309033 - 0.0640878700088575im, 59.999687034378844 - 0.0806813915285931im, 59.9995039874537 - 0.10157094300450688im, 59.99921388448605 - 0.12786842350801764im, 59.998754119492965 - 0.1609731809876083im, 59.99802548435043 - 0.20264599230939304im, 59.99687078413768 - 0.25510175325029183im, 59.99504098053324 - 0.32112528031798987im, 59.992141622554456 - 0.40421532690447065im, 59.9875481704185 - 0.5087623642730564im, 59.9802723581268 - 0.640265382273616im, 59.968751808224795 - 0.8055908977318887im, 59.95052013417703 - 1.0132714385648423im, 59.921692918682695 - 1.2738268988855665im, 59.8761749792913 - 1.60006250698322im, 59.80445809740748 - 2.007237410438142im, 59.6918477293066 - 2.5128837713225196im, 59.51596999235755 - 3.1358521210035586im, 59.24356289157821 - 3.8938296223635462im, 58.82705380048551 - 4.798148020120478im, 58.2026072157431 - 5.844424496272275im, 57.2934629935298 - 6.998304150627912im, 56.024771906796715 - 8.179007460670997im, 54.3552245003153 - 9.250950107344007im, 52.31992912685222 - 10.041559893562525im, 50.05773309668558 - 10.39701137310214im, 47.78711585080809 - 10.254876828903228im, 45.72874597195891 - 9.680443485225494im, 44.02610227843614 - 8.834433869633957im, 42.717590780981574 - 7.9013403159944975im, 41.76134701515096 - 7.032205209856159im, 41.07739917834309 - 6.325170132317829im, 40.5790713626408 - 5.833136293594482im, 40.18691233567748 - 5.579991566176524im, 39.83001443928053 - 5.5741873929073im, 39.44102848047991 - 5.815888181231117im, 38.94935222500849 - 6.2971344897657255im, 38.275974605127004 - 6.995017874645803im, 37.334371044106035 - 7.858360737189663im, 36.043696324073395 - 8.791435267919733im, 34.35937973204408 - 9.645520261937902im, 32.315240491466035 - 10.236685898139323im, 30.049973919645815 - 10.401430139865122im, 27.782154253332997 - 10.06878255632554im, 25.732360316897502 - 9.295516623144646im, 24.044192675661563 - 8.232748343815807im, 22.757273385145087 - 7.053695006993734im, 21.832890680497236 - 5.89634542049547im, 21.196921596987035 - 4.843979180345442im, 20.77224780910884 - 3.932761972354752im, 20.49429310095126 - 3.1681232549055567im, 20.314746245860412 - 2.5392190279952804im, 20.199750615608323 - 2.0285165834724554im, 20.126499999058048 - 1.6171486094290157im, 20.080002695308444 - 1.2874917366376704im, 20.050552855488363 - 1.0241726410295686im, 20.031926452850442 - 0.8142735977388792im, 20.020156094585083 - 0.6471741372561152im, 20.012722363867915 - 0.5142561274993461im, 20.00802915229575 - 0.4085821602267157im, 20.00506680268277 - 0.32459547484470364im]
  f = [0.001, 0.0012589254117941675, 0.001584893192461114, 0.001995262314968879, 0.0025118864315095794, 0.0031622776601683794, 0.003981071705534973, 0.005011872336272725, 0.00630957344480193, 0.007943282347242814, 0.01, 0.012589254117941675, 0.015848931924611134, 0.0199526231496888, 0.025118864315095794, 0.03162277660168379, 0.039810717055349734, 0.05011872336272723, 0.06309573444801933, 0.07943282347242814, 0.1, 0.12589254117941673, 0.15848931924611132, 0.19952623149688797, 0.251188643150958, 0.31622776601683794, 0.3981071705534972, 0.5011872336272722, 0.6309573444801932, 0.7943282347242815, 1.0, 1.2589254117941673, 1.5848931924611136, 1.9952623149688795, 2.51188643150958, 3.1622776601683795, 3.9810717055349722, 5.011872336272722, 6.309573444801933, 7.943282347242816, 10.0, 12.589254117941675, 15.848931924611133, 19.952623149688797, 25.118864315095795, 31.622776601683793, 39.810717055349734, 50.11872336272722, 63.09573444801933, 79.43282347242814, 100.0, 125.89254117941675, 158.48931924611142, 199.52623149688787, 251.18864315095797, 316.2277660168379, 398.1071705534973, 501.18723362727246, 630.957344480193, 794.3282347242814, 1000.0]
  return f, Z
end

function extrapol_almost_peak()
  # case with almost semi-semi-circle
  Z = Any[59.9998025309033 - 0.0640878700088575im, 59.999687034378844 - 0.0806813915285931im, 59.9995039874537 - 0.10157094300450688im, 59.99921388448605 - 0.12786842350801764im, 59.998754119492965 - 0.1609731809876083im, 59.99802548435043 - 0.20264599230939304im, 59.99687078413768 - 0.25510175325029183im, 59.99504098053324 - 0.32112528031798987im, 59.992141622554456 - 0.40421532690447065im, 59.9875481704185 - 0.5087623642730564im, 59.9802723581268 - 0.640265382273616im, 59.968751808224795 - 0.8055908977318887im, 59.95052013417703 - 1.0132714385648423im, 59.921692918682695 - 1.2738268988855665im, 59.8761749792913 - 1.60006250698322im, 59.80445809740748 - 2.007237410438142im, 59.6918477293066 - 2.5128837713225196im, 59.51596999235755 - 3.1358521210035586im, 59.24356289157821 - 3.8938296223635462im, 58.82705380048551 - 4.798148020120478im, 58.2026072157431 - 5.844424496272275im, 57.2934629935298 - 6.998304150627912im, 56.024771906796715 - 8.179007460670997im, 54.3552245003153 - 9.250950107344007im, 52.31992912685222 - 10.041559893562525im, 50.05773309668558 - 10.39701137310214im, 47.78711585080809 - 10.254876828903228im, 45.72874597195891 - 9.680443485225494im, 44.02610227843614 - 8.834433869633957im, 42.717590780981574 - 7.9013403159944975im, 41.76134701515096 - 7.032205209856159im, 41.07739917834309 - 6.325170132317829im, 40.5790713626408 - 5.833136293594482im, 40.18691233567748 - 5.579991566176524im, 39.83001443928053 - 5.5741873929073im, 39.44102848047991 - 5.815888181231117im, 38.94935222500849 - 6.2971344897657255im, 38.275974605127004 - 6.995017874645803im, 37.334371044106035 - 7.858360737189663im, 36.043696324073395 - 8.791435267919733im, 34.35937973204408 - 9.645520261937902im, 32.315240491466035 - 10.236685898139323im]
  f = [0.001, 0.0012589254117941675, 0.001584893192461114, 0.001995262314968879, 0.0025118864315095794, 0.0031622776601683794, 0.003981071705534973, 0.005011872336272725, 0.00630957344480193, 0.007943282347242814, 0.01, 0.012589254117941675, 0.015848931924611134, 0.0199526231496888, 0.025118864315095794, 0.03162277660168379, 0.039810717055349734, 0.05011872336272723, 0.06309573444801933, 0.07943282347242814, 0.1, 0.12589254117941673, 0.15848931924611132, 0.19952623149688797, 0.251188643150958, 0.31622776601683794, 0.3981071705534972, 0.5011872336272722, 0.6309573444801932, 0.7943282347242815, 1.0, 1.2589254117941673, 1.5848931924611136, 1.9952623149688795, 2.51188643150958, 3.1622776601683795, 3.9810717055349722, 5.011872336272722, 6.309573444801933, 7.943282347242816, 10.0, 12.589254117941675]
  return f, Z
end

function extrapol_hard_case()
  f = [0.001, 0.0012589254117941675, 0.001584893192461114, 0.001995262314968879, 0.0025118864315095794, 0.0031622776601683794, 0.003981071705534973, 0.005011872336272725, 0.00630957344480193, 0.007943282347242814, 0.01, 0.012589254117941675, 0.015848931924611134, 0.0199526231496888, 0.025118864315095794, 0.03162277660168379, 0.039810717055349734, 0.05011872336272723, 0.06309573444801933, 0.07943282347242814, 0.1, 0.12589254117941673, 0.15848931924611132, 0.19952623149688797, 0.251188643150958, 0.31622776601683794, 0.3981071705534972, 0.5011872336272722, 0.6309573444801932, 0.7943282347242815, 1.0, 1.2589254117941673, 1.5848931924611136, 1.9952623149688795, 2.51188643150958, 3.1622776601683795]
  Z = Any[59.9998025309033 - 0.0640878700088575im, 59.999687034378844 - 0.0806813915285931im, 59.9995039874537 - 0.10157094300450688im, 59.99921388448605 - 0.12786842350801764im, 59.998754119492965 - 0.1609731809876083im, 59.99802548435043 - 0.20264599230939304im, 59.99687078413768 - 0.25510175325029183im, 59.99504098053324 - 0.32112528031798987im, 59.992141622554456 - 0.40421532690447065im, 59.9875481704185 - 0.5087623642730564im, 59.9802723581268 - 0.640265382273616im, 59.968751808224795 - 0.8055908977318887im, 59.95052013417703 - 1.0132714385648423im, 59.921692918682695 - 1.2738268988855665im, 59.8761749792913 - 1.60006250698322im, 59.80445809740748 - 2.007237410438142im, 59.6918477293066 - 2.5128837713225196im, 59.51596999235755 - 3.1358521210035586im, 59.24356289157821 - 3.8938296223635462im, 58.82705380048551 - 4.798148020120478im, 58.2026072157431 - 5.844424496272275im, 57.2934629935298 - 6.998304150627912im, 56.024771906796715 - 8.179007460670997im, 54.3552245003153 - 9.250950107344007im, 52.31992912685222 - 10.041559893562525im, 50.05773309668558 - 10.39701137310214im, 47.78711585080809 - 10.254876828903228im, 45.72874597195891 - 9.680443485225494im, 44.02610227843614 - 8.834433869633957im, 42.717590780981574 - 7.9013403159944975im, 41.76134701515096 - 7.032205209856159im, 41.07739917834309 - 6.325170132317829im, 40.5790713626408 - 5.833136293594482im, 40.18691233567748 - 5.579991566176524im, 39.83001443928053 - 5.5741873929073im, 39.44102848047991 - 5.815888181231117im]
  return f, Z
end