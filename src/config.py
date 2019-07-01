import os

src = os.path.dirname(__file__)
ROOT = src.rsplit(os.sep, 1)[0]

##############################################################################
# Common File Directories
# input_dir = os.path.join(root, 'data', 'input')
input_dir = os.path.join(ROOT, 'data', 'input')
pdb_files_dir = os.path.join(input_dir, 'pdb_files')

store_dir = os.path.join(ROOT, 'data', 'store')
output_dir = os.path.join(ROOT, 'data', 'output')

##############################################################################

# Parameters
offsets = (-8, 22)
HC_THRESHOLD = 5.
CC_THRESHOLD = 2.8

_aa_index = [('ALA', 'A'),
             ('CYS', 'C'),
             ('ASP', 'D'),
             ('GLU', 'E'),
             ('PHE', 'F'),
             ('GLY', 'G'),
             ('HIS', 'H'),
             ('HSE', 'H'),
             ('HSD', 'H'),
             ('ILE', 'I'),
             ('LYS', 'K'),
             ('LEU', 'L'),
             ('MET', 'M'),
             ('MSE', 'M'),
             ('ASN', 'N'),
             ('PRO', 'P'),
             ('GLN', 'Q'),
             ('ARG', 'R'),
             ('SER', 'S'),
             ('THR', 'T'),
             ('VAL', 'V'),
             ('TRP', 'W'),
             ('TYR', 'Y')]

AA3_to_AA1 = dict(_aa_index)

prosite_extract_path = os.path.join(input_dir, "prosite_extract.txt")
ioncom_path = os.path.join(input_dir, "allid_reso3.0_len50_nr40.txt")

pname_cid_path = os.path.join(store_dir, "pname_cid_map.pkl")

pdb_folder = os.path.join(store_dir, "pdb_files")
pdb_folder_ioncom = os.path.join(store_dir, "pdb_files_ioncom")
fasta_fpath = os.path.join(store_dir, "prosite_seqs.fasta")

meme_out = os.path.join(store_dir, "meme_out")
meme_txt_path = os.path.join(meme_out, "meme.txt")

mast_out = os.path.join(store_dir, "mast_out")
mast_txt_path = os.path.join(mast_out, "mast.txt")

ptr_data_path = os.path.join(store_dir, "ptr_data.pkl")

prosite_pdb_list = [
    '1A03', '1A29', '1A2X', '1A75', '1AHR', '1AJ4', '1AJ5', '1AK8', '1ALV',
    '1ALW', '1AP4', '1AUI', '1AVS', '1B1G', '1B4C', '1B7T', '1B8C', '1B8L',
    '1B8R', '1B9A', '1BJF', '1BLQ', '1BMO', '1BOC',  '1BOD', '1BU3', '1C07',
    '1C7V', '1C7W', '1CB1', '1CDL', '1CDM', '1CDN', '1CDP', '1CFC', '1CFD',
    '1CFF', '1CFP', '1CKK', '1CLB', '1CLL', '1CLM', '1CM1', '1CM4', '1CMF',
    '1CMG', '1CNP', '1CTA',  '1CTD', '1CTR', '1DEG', '1DF0', '1DFK', '1DFL',
    '1DGU', '1DGV',  '1DJG', '1DJH', '1DJI', '1DJW', '1DJX', '1DJY', '1DJZ',
    '1DMO', '1DT7', '1DTL', '1DVI', '1EH2', '1EJ3', '1EL4', '1EXR', '1F4O',
            '1F4Q', '1F54', '1F55', '1F70', '1F71', '1F8H', '1FF1', '1FI5',
            '1FI6', '1FPW', '1FW4', '1G33', '1G4Y', '1G8I', '1GGW', '1GGZ',
            '1GJY', '1H4B', '1HQV', '1HT9', '1I84', '1IG5', '1IGV', '1IH0',
            '1IJ5', '1IJ6', '1IKU', '1IQ3', '1IQ5', '1IRJ', '1IWQ', '1J1D',
            '1J1E', '1J7O', '1J7P', '1JBA', '1JC2', '1JF0', '1JF2', '1JFJ',
            '1JFK', '1JSA', '1JUO', '1JWD', '1K2H', '1K8U', '1K90', '1K93',
            '1K94', '1K95', '1K96', '1K9K', '1K9P', '1K9U', '1KCY', '1KFU',
            '1KFX', '1KK7', '1KK8', '1KQM', '1KQV', '1KSM', '1KWO', '1L2O',
            '1L7Z', '1LA0', '1LA3', '1LIN', '1LKJ', '1LVC', '1LXF', '1M31',
            '1M39', '1M63', '1M8Q', '1MF8', '1MHO', '1MQ1', '1MR8', '1MUX',
            '1MVW', '1MWN', '1MXE', '1MXL', '1N0Y', '1N65', '1NCX', '1NCY',
            '1NCZ', '1NIW', '1NP8', '1NPQ', '1NSH', '1NUB', '1NWD', '1NX0',
            '1NX1', '1NX2', '1NX3', '1NYA', '1O18', '1O19', '1O1A', '1O1B',
            '1O1C', '1O1D', '1O1E', '1O1F', '1O1G', '1OHZ', '1OMD', '1OMR',
            '1OMV', '1OOJ', '1OQP', '1OSA', '1OZS', '1PAL', '1PK0', '1PON',
            '1PRW', '1PSB', '1PSR', '1PVA', '1PVB', '1Q80', '1QIV', '1QIW',
            '1QLK', '1QLS', '1QV0', '1QV1', '1QVI', '1QX2', '1QX5', '1QX7',
            '1QXP', '1REC', '1RFJ', '1RJV', '1RK9', '1RRO', '1RTP', '1RWY',
            '1S1E', '1S26', '1S36', '1S3P', '1S5G', '1S6C', '1S6I', '1S6J',
            '1SBJ', '1SCM', '1SCV', '1SK6', '1SKT', '1SL7', '1SL8', '1SL9',
            '1SMG', '1SNL', '1SPY', '1SR6', '1SRA', '1SW8', '1SY9', '1SYM',
            '1TCF', '1TCO', '1TIZ', '1TN4', '1TNP', '1TNQ', '1TNW', '1TNX',
            '1TOP', '1TRF', '1TTX', '1U5I', '1UHH', '1UHI', '1UHJ', '1UHK',
            '1UP5', '1UWO', '1WDC', '1WRK', '1WRL', '1WRZ', '1X02', '1XA5',
            '1XFU', '1XFV', '1XFW', '1XFX', '1XFY', '1XFZ', '1XK4', '1XO5',
            '1XVJ', '1XYD', '1Y0V', '1Y1A', '1Y6W', '1YR5', '1YRT', '1YRU',
            '1YTZ', '1YV0', '1YX7', '1YX8', '1ZAC', '1ZFS', '1ZMZ', '1ZOT',
            '1ZUZ', '2A4J', '2AAO', '2AMI', '2B1U', '2B59', '2BBM', '2BBN',
            '2BCA', '2BCB', '2BCX', '2BE4', '2BE6', '2BEC', '2BKH', '2BKI',
            '2BL0', '2CCL', '2CNP', '2COL', '2CT9', '2CTN', '2D8N', '2DFS',
            '2DOQ', '2E6W', '2F2O', '2F2P', '2F33', '2F3Y', '2F3Z', '2F8P',
            '2FOT', '2G9B', '2GGM', '2GGZ', '2GV5', '2H61', '2HET', '2HF5',
            '2HQW', '2I08', '2I18', '2I2R', '2I94', '2ISD', '2IX7', '2JC2',
            '2JPT', '2JQ6', '2JT0', '2JT3', '2JT8', '2JTT', '2JTZ', '2JU0',
            '2JUL', '2JWW', '2JXC', '2JXL', '2JZI', '2K0E', '2K0F', '2K0J',
            '2K2F', '2K2I', '2K3S', '2K61', '2K7B', '2K7C', '2K7D', '2K7O',
            '2KAX', '2KAY', '2KBM', '2KDH', '2KDU', '2KFF', '2KFG', '2KFH',
            '2KFX', '2KGB', '2KGR', '2KHN', '2KNE', '2KQY', '2KRD', '2KSP',
            '2KUG', '2KUH', '2KXW', '2KYC', '2KYF', '2KZ2', '2L0P', '2L1R',
            '2L2E', '2L4H', '2L4I', '2L50', '2L51', '2L53', '2L7L', '2L98',
            '2LAN', '2LAP', '2LCP', '2LGF', '2LHH', '2LHI', '2LHL', '2LL6',
            '2LL7', '2LLO', '2LLQ', '2LLS', '2LLT', '2LLU', '2LM5', '2LMT',
            '2LMU', '2LMV', '2LNK', '2LP2', '2LP3', '2LQC', '2LQP', '2LUC',
            '2LUX', '2LV6', '2LV7', '2LVI', '2LVJ', '2LVK', '2LVV', '2M0J',
            '2M0K', '2M1K', '2M28', '2M29', '2M3S', '2M3W', '2M49', '2M55',
            '2M5E', '2M7K', '2M7M', '2M7N', '2MA2', '2MAZ', '2MBX', '2MES',
            '2MG5', '2MGU', '2MKP', '2MLE', '2MLF', '2MRD', '2MYS', '2MZP',
            '2N27', '2N6A', '2N77', '2N79', '2N7L', '2N8J', '2N8Y', '2N8Z',
            '2NA0', '2NCO', '2NCP', '2NLN', '2NXQ', '2NZ0', '2O5G', '2O60',
            '2OBH', '2OPO', '2P6B', '2PAL', '2PAS', '2PMY', '2PQ3', '2PRU',
            '2PSR', '2PVB', '2Q4U', '2Q91', '2QPT', '2R28', '2R2I', '2RGI',
            '2RO9', '2RRT', '2SAS', '2SCP', '2TN4', '2V01', '2V02', '2V53',
            '2VAS', '2VAY', '2VB6', '2VN5', '2VN6', '2VRG', '2W49', '2W4A',
            '2W4G', '2W4H', '2W4T', '2W4U', '2W4V', '2W4W', '2W73', '2WEL',
            '2WND', '2WOR', '2WOS', '2X0G', '2X51', '2Y4V', '2YGG', '2ZN8',
            '2ZN9', '2ZND', '2ZNE', '2ZRS', '2ZRT', '3A4U', '3A8R', '3AAJ',
            '3AAK', '3B32', '3BOW', '3BXK', '3BXL', '3BYA', '3C1V', '3CGA',
            '3CLN', '3CR2', '3CR4', '3CR5', '3CS1', '3CTN', '3CZT', '3D0Y',
            '3D10', '3DD4', '3DF0', '3DVE', '3DVJ', '3DVK', '3DVM', '3E3R',
            '3EK4', '3EK7', '3EK8', '3EKH', '3EVR', '3EVU', '3EVV', '3EWT',
            '3EWV', '3F45', '3FS7', '3FWB', '3FWC', '3G43', '3GK1', '3GK2',
            '3GK4', '3GN4', '3GOF', '3GP2', '3H4S', '3HCM', '3HR4', '3I5F',
            '3I5G', '3I5H', '3I5I', '3ICB', '3IF7', '3IFK', '3IQO', '3IQQ',
            '3J04', '3J41', '3JTD', '3JVT', '3K21', '3KCP', '3KF9', '3KO0',
            '3L9I', '3LCP', '3LI6', '3LK0', '3LK1', '3LL8', '3LLE', '3M0W',
            '3NXA', '3O77', '3O78', '3OX5', '3OX6', '3OXQ', '3PAL', '3PAT',
            '3PM8', '3PSR', '3PX1', '3QJK', '3QRX', '3RLZ', '3RM1', '3RV5',
            '3SG2', '3SG3', '3SG4', '3SG5', '3SG6', '3SG7', '3SJQ', '3SUI',
            '3UCT', '3UCW', '3UCY', '3ULG', '3WFN', '3WHT', '3WHU', '3WLC',
            '3WLD', '3WNX', '3WXA', '3ZWH', '4ANJ', '4AQI', '4AQJ', '4AQR',
            '4BW7', '4BW8', '4BYA', '4BYF', '4C0J', '4C0K', '4C0L', '4CFQ',
            '4CFR', '4CID', '4CLN', '4CPV', '4DBP', '4DBQ', '4DCK', '4DIR',
            '4DJC', '4DS7', '4DUQ', '4E50', '4E53', '4EHQ', '4ETO', '4F0Z',
            '4FL4', '4FQO', '4G27', '4G28', '4GGF', '4GOW', '4GUK', '4HEX',
            '4HSZ', '4I2Y', '4I5J', '4I5K', '4I5L', '4I5N', '4ICB', '4IL1',
            '4J9Y', '4J9Z', '4JPZ', '4JQ0', '4L79', '4L9M', '4LZX', '4M1L',
            '4M2O', '4M2P', '4M2Q', '4MBE', '4MEW', '4MLW', '4MRX', '4MRY',
            '4MSP', '4MVF', '4N1F', '4N1G', '4N5X', '4NQG', '4NSC', '4NSD',
            '4OKH', '4OR9', '4ORA', '4ORB', '4ORC', '4OV2', '4OVN', '4OY4',
            '4P2Y', '4P5X', '4P60', '4PAL', '4PCW', '4PDZ', '4PE0', '4PE1',
            '4PE4', '4PE7', '4PHJ', '4PHK', '4PHM', '4PHN', '4PJJ', '4Q57',
            '4Q5U', '4QNH', '4QOX', '4R8G', '4RGJ', '4RJD', '4TNC', '4U8D',
            '4UMO', '4UPG', '4UPU', '4USL', '4V0C', '4WQ2', '4WQ3', '4XJK',
            '4XYN', '4Y99', '4YBH', '4YGB', '4YGC', '4YGD', '4YGE', '4YI8',
            '4YI9', '4YRU', '4ZCU', '4ZCV', '4ZLK', '5A2H', '5AEQ', '5AER',
            '5AFP', '5COC', '5CPV', '5CSF', '5CSI', '5CSJ', '5CSN', '5D43',
            '5D69', '5D7F', '5DBR', '5DKN', '5DKQ', '5DKR', '5DOW', '5DSU',
            '5E1K', '5E1N', '5E1P', '5ER4', '5ER5', '5G4P', '5G58', '5G5D',
            '5GGM', '5GQQ', '5H53', '5H7D', '5HIT', '5HLO', '5HLV', '5HYD',
            '5I0I', '5I8N', '5J03', '5J7J', '5J8H', '5JJG', '5JOJ', '5JOL',
            '5JQA', '5JTH', '5K7L', '5K89', '5K8Q', '5KSZ', '5KTY', '5KU1',
            '5LPU', '5M6C', '5MRA', '5NIN', '5OEO', '5OTJ', '5PAL', '5SVE',
            '5T0X', '5T7C', '5TBY', '5TNC', '5TP6', '5V02', '5V03', '5V7X',
            '5VE9', '5VLN', '5VMS', '5W1F', '5W88', '5WBX', '5WC5', '5WCL',
            '5WSU', '5WSV', '5XND', '5ZAB', '5ZGM', '6AGI', '6AGJ', '6ALE',
            '6B8L', '6B8M', '6B8N', '6B8P', '6B8Q', '6BNV', '6C1D', '6C1G',
            '6C1H', '6CNM', '6CNN', '6CNO', '6CZQ', '6DAD', '6DAE', '6DAF',
            '6DAH', '6DMW', '6DS2', '6E2F', '6E2G', '6FEG', '6FEH', '6FIE',
            '6GDK', '6GDL', '6HCS', '6IIE', '6MV3']