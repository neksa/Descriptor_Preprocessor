import urllib.parse
import urllib.request


def convert(from_id, to_id, input_ids):
    """
    Assumption that new_id is unique, and that for multiple matches to
    old_id, the first is selected. Works well when new_id is PDB_ID.
    :param from_id:
    :param to_id: "PDB_ID"
    :param input_ids:
    :return:
        # ['MSHA_SALTO\tmshA', 'NADE_STRCO\tnadE']
    """
    assert isinstance(input_ids, list)
    query_line = " ".join(input_ids)
    url = 'https://www.uniprot.org/uploadlists/'

    params = {'from': from_id, 'to': to_id, 'format': 'tab',
              'query': query_line}
    # params = {'from': 'ACC+ID', 'to': 'PDB_ID', 'format': 'tab',
    #           'query': 'P30575 P04764 P07323 P17182 P17183 P15429 P21550'}
    # params = {'from': 'ACC', 'to': 'PDB_ID', 'format': 'tab',
    #           'query': 'P30575 P04764 P07323 P17182 P17183 P15429 P21550 '
    #                    'P06733 Q9C9C4 Q5XD01 Q97QS2 P13929 P00925 P25696 '
    #                    'P09104 P00924 P0A6P9 Q8IJN7 O51312 Q8DTS9 P42848 '
    #                    'Q9XSJ4 Q2G028 P25704 Q27727 Q473G4 P69951 B2IPX8 '
    #                    'Q1CV03 Q4FM37 B1LCF4 B0SRL5 B0K882 B1ZEJ8 O52191 '
    #                    'B4SGB2 Q7NAY0 A1AEW7 Q9UAL5 P64080 P42040 C1CXJ3 '
    #                    'Q02RA7 C1AY93 B3Q6L0 A9VQ48 C4L5H4 A9HJ75 Q9PJF3 '
    #                    'A8H1S5 A6UUM2 B2TPW4 Q9Z7A6 Q9A7J9 B9KEN0 A5IQY0 '
    #                    'B1IDB9 Q3JCT1 A7G9Y3 A4SRC1 B2A6Z1 Q5FH95 A1K7F6 '
    #                    'A3N1B9 B7KRB4 Q0BSX3 A5FN12 O29133 B5EGF2 A2SAU7 '
    #                    'A8G9W1 Q9HDT3 A5FRM5 Q3AFC8 A8LQL4 A6L3M9 Q6F0Z7 '
    #                    'A3PIV3 A0KGH3 Q49W03 B0S1G7 A0K8N1 B6EKL8 B0SYX8 '
    #                    'B5YF30 B5Z3E3 A6TZQ5 Q6GIL4 B7IFN4 Q2SXC5 Q8DL40 '
    #                    'A4W2T1 Q9PDT8 C1CEB3 Q6AAB8 A6SXG3 B8HQY9 P69949 '
    #                    'B1V8R6 Q886M3 B4RBW1 B1LZU6 B2IKR4 Q07ND9 B1KPT6 '
    #                    'A5VQQ7 A0B7E8 B0TK04 Q39EV9 B3EEQ2 B6JFY2 C1A306 '
    #                    'Q9U615 Q7N835 Q8L202 B5Z9S9 Q21V37 B5BF02 Q5LQL4 '
    #                    'O85348 B3EL51 Q2SKX0 Q2JIT3 Q8R967 B3GY00 Q89KV6 '
    #                    'Q89Z05 Q9JZ53 Q12560 B8GW68 Q65EN2 A9KQ46 A4QCV0 '
    #                    'B8E0W1 A0LW71 Q47DI1 Q31G68 C3NE30 B8J467 Q0C0R3 '
    #                    'A9R1D1 P86210 Q6BTB1 B3CNP3 B1JK09 P64074 B9MS47 '
    #                    'Q1I646 B8DDA1 B5ZAZ0 Q1MS78 Q04KG2 Q7V377 B6J3X4 '
    #                    'Q96VP4 Q1CLT2 P42894 Q1JHQ6 Q31CX3 Q9W7L1 Q2IWD0 '
    #                    'B5E4P1 Q1J7I5 Q8ZBN2 Q27527 P0A6Q2 Q1C3Y6 Q87DY6 '
    #                    'B7JW48 Q2JRS2 Q48F79 B4T483 Q6FQY4 B1MKG1 A5V1W0 '
    #                    'A3QC77 A1A143 Q5LG64 A6US08 A1UL09 Q9XDS7 Q8CPY3 '
    #                    'P57492 B1XLD0 B1VFL0 A8FNY3 P40370 Q31XL1 B2U9C3 '
    #                    'C3LDI0 Q0A7K4 Q65VZ7 A4YVC0 P64079 Q8ZYE7 B4RMD8 '
    #                    'Q72F92 B7LEJ8 B4RVU5 A5FUW3 Q8NKC2 A2RIX1 P31683 '
    #                    'Q0B080 Q3IQT0 A7HNS8 Q1WSY0 Q057H3 A1SF66 B9E6B1 '
    #                    'A6LR09 P0CX11 Q1JML5 Q6G173 Q049Y3 B2VFY8 Q2L088 '
    #                    'Q54RK5 Q050L5 B1JB38 Q66ED8 Q6ADR6 B9LPW6 B0BA40 '
    #                    'A4WR18 C3MPU0 A5U166 Q3K2B2 Q136E5 Q215A2 Q72H85 '
    #                    'Q870B9 Q2FTN3 O26149 A5IYA8 B2TZF4 Q8TUV6 P47647 '
    #                    'Q9LEI9 B2SXU9 A3PAS6 C3PFA7 B1LQB2 Q7W5N9 Q4R5L2 '
    #                    'Q8Y0B5 Q6HBF3 A9WCM4 A0JU21 P42448 A1S4D7 B2GM13 '
    #                    'Q03SL5 C4ZZT2 A8IBF3 B7JFG3 A8AY46 Q8UFH1 Q5KVE7 '
    #                    'Q821H7 B0RH60 C3MBJ9 B0BQ53 B4U2B8 Q6LMT1 Q5L4S5 '
    #                    'B8CJP6 O02654 Q2NVN7 A0M568 C1C7C0 Q9LEJ0 P08734 '
    #                    'A4F804 Q31QJ8 Q9HJT1 Q97ZJ3 B8HEU5 Q74K78 Q83B44 '
    #                    'Q83H73 Q5X3L4 Q1KYT0 Q3JQQ6 Q47WR1 Q042F4 B0RU05 '
    #                    'B4TTY5 A9BIS7 Q04SI9 Q03LI0 A0RMH2 C0Z6L3 Q1J2H6 '
    #                    'A4SXE9 A6GZ69 B0KSB9 O84591 B1MVW3 A1VLH1 A6VUU9 '
    #                    'A1APJ8 C5B8X2 Q0BKV1 Q8XKU4 A3NB82 A0QBX4 Q8G0G3 '
    #                    'C5C093 A1JJR4 A9IS50 P42897 Q7RA60 B1IU62 Q4KHF6 '
    #                    'Q3KLB0 A7INB6 A4YHC1 Q3J3H9 Q1LTN8 B7N715 Q5N3P4 '
    #                    'Q5FNN5 A4Y943 A1R485 A5GPE0 P9WNL0 Q98Q50 A1B9D2 '
    #                    'C1FQW6 B0B8G1 B7HW52 Q5HSC1 Q4QLX6 Q8ENP5 Q5HB46 '
    #                    'P42896 Q661T0 Q88YH3 A5V3E8 B3EUH8 Q0KCE2 E6ER18 '
    #                    'Q9W7L2 Q8PLS0 Q5GTG4 B1IBR3 P07322 B2SUA6 Q9PQV9 '
    #                    'Q48UF7 A4JFY5 B2K561 P19140 Q4J920 Q2IID9 A2BP03 '
    #                    'B1XUR8 Q7VDY0 Q6RG04 B0S8S8 A8G2L5 O74286 P64077 '
    #                    'B7J1R2 Q3AVW5 Q8RP81 B0UV89 A5IIQ6 Q2FLB5 Q7V483 '
    #                    'C5A2S7 Q74IV0 B2S455 C1AM14 B1YT09 Q7NIR1 Q3YRX9 '
    #                    'P0A6Q0 Q2A278 Q8EBR0 Q9CHS7 A0L7X8 Q88MF9 B7V7V4 '
    #                    'Q12PZ4 P64081 B5FTU9 P15007 A6VJ31 A7H622 Q0VQD6 '
    #                    'Q21LC2 B4SR86 B2I2A5 O69174 A5N2N5 Q2FIL7 A5EK13 '
    #                    'A7HXW7 A8MFY2 B8D7U7 P59566 A6U8E5 Q5NZ69 B7GTK2 '
    #                    'C5D7M1 P26301 B3E2S8 Q9HQI9 P0CX10 Q2LR33 P51555 '
    #                    'Q9HXZ5 B4U9X7 O57391 B1AIH2 Q8D2K1 Q1JCN8 Q87WD5 '
    #                    'Q73HQ2 Q4ZWE0 Q898R0 Q9F3P9 Q24YW4 A1RHF3 Q1B439 '
    #                    'Q0BDS3 Q8KB35 B5QW40 Q8EW32 Q76KF9 C6DDJ5 Q7RV85 '
    #                    'B2S5Y3 Q60173 B9KQU4 A7I455 Q2W698 A5WG13 Q0I1Z1 '
    #                    'Q1QMI9 Q57D07 Q181T5 B0JMP6 Q3M7B2 P26300 A4WDW7 '
    #                    'B9M3M1 B7NV69 A0R3B8 Q0HL72 Q3SRK4 Q74AR6 Q2FQL9 '
    #                    'A1KUB6 Q7VNM6 B5XV19 B2J1R2 B7I918 Q27655 B8F8L5 '
    #                    'B3PJB3 A6TU30 A0RKS3 Q1QZX7 B7H227 Q3ATQ5 B7MLA0 '
    #                    'C6E471 C3KZ49 A8FHJ0 B2URY4 P99088 Q8GE63 Q11HU8 '
    #                    'Q39T27 B8I4U1 P42222 Q8YRB0 Q7VQH3 A9IIP8 A7GUR7 '
    #                    'P0DM31 Q967Y8 C1F9E6 Q7MHQ1 B6YQN5 B9L7U3 Q71WX1 '
    #                    'Q631M2 A1KHG1 B6IQ30 A9AGW2 Q468E2 A0KU82 Q1H011 '
    #                    'Q72QZ8 Q6N5U6 Q032H8 Q0AD93 Q98MZ3 B8GQ75 Q2J619 '
    #                    'C3N5G6 Q3IDM2 Q2YPV0 Q6BI20 A6V1F3 P64076 A2SJR2 '
    #                    'C0QRV6 Q83FF7 Q7UIR2 A4X3B7 B4EDA1 Q7U0U6 A8L149 '
    #                    'A2SSV1 B2S044 A1AW20 Q164E3 A6VR00 Q97L52 B2UY16 '
    #                    'Q9CD42 A1WL86 A5F5I3 Q5B135 Q5NGW8 Q9JU46 Q03AK4 '
    #                    'Q2GGS6 Q9UAE6 A6WR28 C0M6K5 Q8RI55 Q1D401 Q1GVS8 '
    #                    'Q8G5I9 P77972 O32513 A5CX71 B9JEY7 P0DA95 Q972B6 '
    #                    'Q67SV9 I0J1J1 Q88VW2 Q5R143 Q8DPS0 C0R5N4 B8ZPW9 '
    #                    'Q4UTP2 A3DBQ5 C1CKJ0 Q316Q0 Q7M8Q0 Q9PVK2 Q32CD6 '
    #                    'P33675 Q11QE1 A2RFE3 Q5GYK4 B6J4X4 B7KL24 A7MQZ0 '
    #                    'B5RDS5 Q28RE8 Q8P9Z3 Q3KH92 Q73P50 Q2GK24 Q256G5 '
    #                    'Q6KZN3 Q7MTV8 Q0TQZ2 Q6M075 Q2NG02 B6YUB8 A8YUV4 '
    #                    'P56252 Q30P06 P48285 B9EAH1 Q9F2Q3 A1TZ48 B9DRR9 '
    #                    'B1GYL7 B5ZNA1 Q57KH0 Q0TE80 P51913 P64078 A0Q5J9 '
    #                    'B5RRF1 Q6KIB0 A5UI73 A4T6L5 Q0SNH5 A1W1S4 A4SGL6 '
    #                    'B8GAX5 Q8PT81 B0V677 Q2G662 Q9UXZ0 Q5F8Z2 A8EWY0 '
    #                    'A3D795 B5EPM7 Q2Y9P0 Q6FAT9 B3QXY4 A6TD53 A8FSS8 '
    #                    'A1TLS5 Q5YQ30 P0A6Q1 Q82HH5 Q96X30 O66778 A1VGV5 '
    #                    'B7J6R4 Q5M0M5 C4XLR9 Q3ZC09 B5XKM7 B3R499 C6BSL8 '
    #                    'C0MH89 A7FLZ5 Q043Z5 Q9BPL7 A3CMA7 Q3ZX11 A0LEC9 '
    #                    'B8JAC4 C3NHN3 A8ZSY5 A0ALD9 Q2RLT8 Q2K8W9 Q6C1F3 '
    #                    'P42895 A9KDH7 A8M2Y3 Q2GD37 A0PYP4 Q1G9S9 B4EUF7 '
    #                    'Q5PAS6 C5CSV6 Q6W3C0 Q03GW5 Q6D182 Q42971 Q6MPQ2 '
    #                    'Q1R7R4 B7LWP5 Q4A740 Q0RCG8 C4LHL7 A5UDD6 Q96X46 '
    #                    'A9M5E5 Q110V4 Q3A578 Q5M561 Q6YQT9 Q5SME1 Q8U477 '
    #                    'A0PW55 Q493N5 B1JUY6 A1US94 A6X0L8 A9LZL4 B6I6H5 '
    #                    'Q6MTZ2 A9FA52 P75189 Q4AA88 A3M5Y1 Q741U7 P33676 '
    #                    'B7IP20 Q27877 B7UHJ5 Q1ISS7 Q43130 Q756H2 Q8K9E0 '
    #                    'Q9K717 B9J4M4 B0R4Y8 A5EW24 C1A7Y3 B3WCW7 B1XDI9 '
    #                    'C3MUV2 A5IDM2 A9BDH2 A3NX13 Q6AM97 P0DA94 Q47SV1 '
    #                    'B1HVS4 Q1MH36 Q38Y18 B2V9N8 A7HCC4 B5DGQ7 Q128E4 '
    #                    'Q604M4 Q03ZK4 Q8GR70 C1DCC3 Q3BUT0 Q55F83 Q74J64 '
    #                    'Q9W7L0 Q92Q98 C4LBR1 B2I937 C5C987 Q0STE0 Q0W8C9 '
    #                    'Q4ZWQ8 Q01YD1 A6W6X2 A5GEW5 Q5HQV0 Q87LQ0 Q9RR60 '
    #                    'A7ZQM2 Q3B1G7 C4Z1M9 A4G737 A2C038 O59605 B1YLD8 '
    #                    'Q9KPC5 Q2RT60 A1W4R1 Q7U3T1 Q4A8B5 Q2YSE8 A7FQP0 '
    #                    'A4ISP4 A2CDE2 B8FKT4 C1CRM6 B5F4N9 Q5ZTX1 Q0T1P7 '
    #                    'A9N9U9 Q46HG5 Q9ZW34 B4TFZ1 Q8YHF0 Q979Z9 B8DTI9 '
    #                    'B2GAM0 Q2NAQ1 Q2SSR3 Q5WV02 C0PXD5 A9N2F4 Q3AGS4 '
    #                    'Q1BHS0 B3QQT2 Q24MW5 C1EZA0 Q5FKM6 Q9Y927 Q5HHP1 '
    #                    'B6JPQ2 Q0HXH0 P29201 A9KYH0 B5FAF9 B9MEQ4 Q0I6K2 '
    #                    'B0C9F0 Q3YY77 C0QI43 A1WWZ2 Q9CIT0 B1Y4K4 Q0S4I1 '
    #                    'B8EJT3 Q5WDK9 Q601S2 A4XKV0 Q6NI61 P74934 Q8KG25 '
    #                    'P43806 Q4FR74 Q8NRS1 Q14IC0 A7ZFN1 A8A3R4 A3Q5F7 '
    #                    'Q1AXJ9 Q15QR6 B0CGT2 A9MF11 Q1Q9K6 B2JIX0 A7Z8Y1 '
    #                    'Q8KNX9 P57975 B4S487 C0RJA3 A4XWS1 Q815K8 B8D9J5 '
    #                    'Q6FTW6 Q62J10 Q3Z8W4 A7WZT2 B3DSV2 Q8TQ79 Q8DC62 '
    #                    'B7LXJ5 A1V5K2 Q12VE5 Q7VW79 A6QF85 B1KTJ8 C3P0A3 '
    #                    'A1BD13 B9KA91 Q6MEY2 A4FZA6 Q086B0 Q4JU51 C4KH32 '
    #                    'Q7WD75 B2HDI5 B2RLL7 Q8F4T8 A4TPY1 Q8FQS7 A7H0M5 '
    #                    'B4UDM8 Q0APT2 Q63SQ0 A2BUI5 C3LQZ0 Q2S4F8 C4K4K1 '
    #                    'Q43321 A6Q5K4 A1TEE4 Q606T2 Q5E326 Q7VIH4 Q6GB54 '
    #                    'A9NF93 Q4L4K7 A1QQ97 A5CTA8 Q18G62 P64075 B0K6X6 '
    #                    'Q04DH2 Q3SL43 Q2P1K8 Q5JEV6 Q9ZMS6 Q030Y9 A1SSQ7 '
    #                    'Q17YV0 A9W6G9 A6LJF0 Q1LPI5 Q13X07 B5RLS2 Q5PEH4 '
    #                    'C1KY94 Q1GI52 B2KBA5 A7NFI9 P9WNL1 B7MZ75 Q2NJ39 '
    #                    'B2FK88 B9DJJ9 B0VQI4 Q81X78 Q70CP7 A8Z1A4 B8E8T1 '
    #                    'Q72XY5 A3ML77 B9JW75 B7HED2'}
    print(params)
    data = urllib.parse.urlencode(params)

    data = data.encode('utf-8')
    req = urllib.request.Request(url, data)
    with urllib.request.urlopen(req) as f:
        response = f.read()
    parsed_response = response.decode('utf-8')
    old_new_id_str = parsed_response.strip().split("\n")
    print(len(parsed_response))
    if len(old_new_id_str) == 1:  # Only header line
        return dict()
    old_new_id_str = old_new_id_str[1:]
    old_new_id_map = dict()
    for line in old_new_id_str:
        old_id, new_id = line.split("\t")
        if old_id not in old_new_id_map:
            old_new_id_map[old_id] = new_id.lower()
            # old_new_id_map[old_id] = new_id
    return old_new_id_map

#
# old_new_id_map = convert(None, None, None)
# print(len(old_new_id_map))