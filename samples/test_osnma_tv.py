"""
test script for Galileo OSNMA based on [2]

[1] Galileo Open Service Navigation Message Authentication (OSNMA)
    Signal-in-Space Interface Control Document (SIS ICD), October, 2023.

[2] Galileo Open Service Navigation Message Authentication (OSNMA)
    Receiver Guidelines Issue 1.3, January, 2024.

@author Rui Hirokawa

"""

import copy
from cssrlib.osnma import osnma, uOSNMA, taginfo
from cssrlib.gnss import prn2sat, uGNSS
from binascii import unhexlify, hexlify

nma = osnma()

# A.6.3 Associated Navigation Data
# sub-frame WN = 1248, TOW = 345660
# TOW = 345631 - 345660
# E1-B I/NAV page (120 bits) from E02
navmsg = [
    '020EB72D6AB6C9DEBDBF3C87EE63C0',
    'BDFA60936A7475EAAAAA7085944100',
    '040E82001A000E5910006D31F00000',
    'A68070D8F793F7AAAAAA68A91DCAC0',
    '06FFFFFFFF0000011248E089E24A80',
    '8C46B2510E7B56AAAAAA7A56FF0BC0',
    '09A23A5555555555552A0110156B40',
    '85F48E4603DD972AAAAA7B8D5E4100',
    '0AADC96FE3DEC7FDC7FFF0FFBDFF00',
    'A6083B13FF8F46EAAAAA723127CAC0',
    '125A95F82358E0521768C7435B86C0',
    '871C9883F9C0976AAAAA5A59088BC0',
    '148BB93DF4BB237D30163105D10840',
    '851E42134280B3EAAAAA6CFB898100',
    '1010D97DB7D45ADCB5A831D98C00C0',
    '9990055251F6BFEAAAAA6F7FCF0AC0',
    '009555555555555555555555538140',
    '918BC4060303E5EAAAAA7BD6FCCBC0',
    '009555555555555555555555538140',
    '918C503FCEAF48EAAAAA63697AC100',
    '010E96441475A496001DACFBAA0500',
    '86FCF304B10B70EAAAAA44E2740AC0',
    '030EBFF0DC080EC363843AC5344440',
    'A39AFFE8F36C45EAAAAA7BD4640BC0',
    '0554A02C22607F7FC009C0A8C6EA80',
    'AAAAB3FEAB6E0EEAAAAA754D4B0100',
    '009555555555555555555555538140',
    '918E47D7C2DBFAAAAAAA4A4B4DCAC0',
    '1010D97DB7D45ADCB5A831D98C00C0',
    '99902DE23AC0002AAAAA5EC4050BC0',
]

# A.3 NMA Header
nma.nma_header = 0b10000010

# A.4 DSM-KROOT
dsm_kr = '2210492204E060610BDF26D77B5BF8C9CBFCF70422081475FD445DF0FFF8CD88299FA4605800207BFEBEAC55024053F30F7C69B35C15E60800AC3B6FE3ED0639952F7B028D86867445961FFE94FB226BFF7006E0C451EE3F8728C177FB5E130DA4B44BBE7EC29522'
dsm_kr_ = unhexlify(dsm_kr)
did = 0
nma.dsm = {did: dsm_kr_}

# A.4.1 DSM-KROOT interpretation
nma.decode_dsm_kroot(did)
result = (hexlify(nma.alp) == b'610bdf26d77b')
print(f"A.4.1 DSM-KROOT interpretation result={result}")

nma.prn_a = 2  # E2
gst_wn = 1248

sat = prn2sat(uGNSS.GAL, nma.prn_a)

# A.4.2 DSM-KROOT verification
pubk_path = 'OSNMA_PublicKey_test.crt'
nma.load_pubkey_pkid(pubk_path, nma.pkid)
result = nma.verify_root_key()
print(f"A.4.2 DSM-KROOT verification result={result}")

# A.5 TESLA Chain Key
# A.5.1 TESLA Chain Key Interpretation
nma.key_c = unhexlify('01D3E3E2667A0A1894D04BDD98CABA17')  # K120
nma.gst_tow = 349170

# A.5.2 TESLA Chain Key Verification
nma.status |= uOSNMA.ROOTKEY_VERIFIED
result = nma.verify_key_chain(nma.key_c, gst_wn, nma.gst_tow)
print(f"A.5.2 TESLA Chain Key Verification for KROOT result={result}")

# tow = 345660 K3
nma.key_p = unhexlify('2E9C651C410CFF23D370497D28AB4B14')
# tow = 345690 K4
nma.key_c = unhexlify('69C00AA7364237A65EBF006AD8DDBC73')
nma.gst_tow = 345690

nma.status |= uOSNMA.KEYCHAIN_VERIFIED
result = nma.verify_key_chain(nma.key_c, gst_wn, nma.gst_tow)
print(f"A.5.2 TESLA Chain Key Verification for K4 result={result}")

# A.6 Tags
# A.6.1 MACK Message Interpretation
nma.gst_tow_p = nma.gst_tow-30
nma.gst_sf_p = nma.set_gst_sf(gst_wn, nma.gst_tow_p)

# MACK Message E02 at WN = 1248, TOW = 345630 seconds
# l1 = '4DA9D1D763DE4FDE4439ED5A180F765C4FFE3D1B0FE7025D4D0A02CF4947DAFF180C0F97FF3ABD2312C42DC3A3CDB117FAADB83B5F0B6FEA88EB0000'
# MACK Message E02 at WN = 1248, TOW = 345660 seconds
l2 = 'E37BC4F858B08F0726A9000C12057BB238C883024FC8FB247323050F91F230FE4D02CF2F1C4D28C1220F2E9C651C410CFF23D370497D28AB4B140000'

# msg1 = unhexlify(l1)
msg2 = unhexlify(l2)

lt = nma.tag_len_t[nma.ts]  # length of MAC
lk = nma.klen_t[nma.ks]  # length of Key
nma.nt = (480-lk)//(lt+16)  # number of tags
bl = (lt+16)//8

nma.tag = bytearray(msg2[0:bl*nma.nt])
result = nma.verify_maclt()
print(f"A.6.2 Tag Sequence Verification result={result}")

nma.key = nma.key_c
result = nma.verify_macseq()
print(f"A.6.4 MACSEQ Verification result={result}")

# A.6.5 Tags Verification
# A.6.5.1 Tag0 Verification
prn = nma.prn_a
for k in range(len(navmsg)//2):
    msg = unhexlify(navmsg[2*k])+unhexlify(navmsg[2*k+1])
    nav, nma_b = nma.load_gal_inav(msg)
    nma.save_gal_inav(nav, prn, nma.gst_tow-60+(k+1)*2)

    nma.hk[prn-1][k] = nma_b[0]              # HK-ROOT message
    nma.mack[prn-1][k*4:k*4+4] = nma_b[1:5]  # MACK message


nma.subfrm[sat] = copy.copy(nma.subfrm_n[sat])

adkd0 = nma.gen_gal_inavmsg(nma.prn_a)
# print('ADKD0=', hexlify(adkd0))

# tag0 = unhexlify('E37BC4F858')  # tow=345660
ctr = 1
i0 = 7*(ctr-1)
tag0 = msg2[i0+0:i0+5]
adkd = 0
cop = msg2[i0+6] & 0xf

prn_d = nma.prn_a
ctr = 1
tag_ = taginfo(nma.gst_sf_p, prn_d, nma.prn_a, adkd, cop, tag0, ctr, adkd0)

m0 = nma.gen_msg(adkd, prn_d, nma.gst_sf_p, ctr, tag_.navmsg)
tag_c = nma.process_mac(m0)[:5]
print(f"A.6.5.1 Tag0 Verification result={tag_c == tag0}")
# result = nma.verify_navmsg(tag_)
# print(f"ADKD0 result={result}")

# A.6.5.2 ADKD4 Verification
ctr = 3
i0 = 7*(ctr-1)
tag4 = msg2[i0+0:i0+5]
adkd = msg2[i0+6] >> 4

adkd4 = nma.gen_gal_utcmsg()
m4 = nma.gen_msg(adkd, prn_d, tag_.gst_sf, ctr, adkd4)
tag_c = nma.process_mac(m4)[:5]
print(f"A.6.5.2 ADKD4 Verification result={tag_c == tag4}")

# A.6.5.3 ADKD12 Verification
ctr = 5
i0 = 7*(ctr-1)
tag12 = msg2[i0+0:i0+5]
adkd = msg2[i0+6] >> 4
cop = msg2[i0+6] & 0xf
prn_d = nma.prn_a
tag_ = taginfo(nma.gst_sf_p, prn_d, nma.prn_a, adkd, cop, tag0, ctr, adkd0)

# tow = 345990 K14
nma.key = unhexlify('ED2B77AAED7B633B58A69551FC388421')

m12 = nma.gen_msg(adkd, prn_d, nma.gst_sf_p, ctr, tag_.navmsg)
tag_c = nma.process_mac(m12)[:5]
print(f"A.6.5.3 ADKD12 Verification result={tag_c == tag12}")

# DSM-PKR
dsm_pkr_ = '717CBE05D9970CFC9E22D0A43A340EF557624453A2E821AADEAC989C405D78BA06956380BAB0D2C939EC6208151040CCFFCF1FB7156178FD1255BA0AECAAA253F7407B6C5DD4DF059FF8789474061301E1C34881DB7A367A913A3674300E21EAB124EF508389B7D446C3E2ECE8D459FBBD3239A794906F5B1F92469C640164FD87120303B2CE64BC207BDD8BC4DF859187FCB686320D63FFA091410FC158FBB77980EAB8884C0D33D6'
nma.root_mt = unhexlify(
    'A10C440F3AA62453526DB4AF76DF8D9410D35D8277397D7053C700D192702B0D')
did = 13
nma.dsm[did] = unhexlify(dsm_pkr_)
result = nma.decode_dsm_pkr(did)
print(f"A.7.2 DSM-PKR Verification result={result}")
