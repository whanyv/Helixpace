from dataclasses import dataclass
import os
import pydicom
import numpy as np
from tqdm import tqdm
import sys
import glob
import matplotlib.pyplot as plt


class primerEncoder:
    def __init__(self):
        pass
        

    def str2bin(str_data):
        res = ''.join(format(ord(i),'08b') for i in str_data)
        return res


    def Bin2dec(binary):
        string = int(binary,2)
        return string


    def bin2str(bin_data):
        str_data = ''
        for i in range(0, len(bin_data), 8):
            temp_data = bin_data[i:i + 8]
            dec_data = primerEncoder.Bin2dec(temp_data)
            str_data = str_data + chr(dec_data)

        return str_data


    def bin2seq(bin_data):
        seq = ''
        for i in range(0, len(bin_data), 2):
            temp = bin_data[i:i+2]
            if temp == '00':
                seq = seq + 'A'
            if temp == '01':
                seq = seq + 'T'
            if temp == '10':
                seq = seq + 'C'
            if temp == '11':
                seq = seq + 'G'
        return seq


    def seq2bin(seq_data):
        bin = ''
        for i in len(seq_data):
            if i == 'A':
                bin = bin + '00'
            if i == 'T':
                bin = bin + '01'
            if i == 'C':
                bin = bin + '10'
            if i == 'G':
                bin = bin + '11'
        return bin


    #actual primer encoding function. primer lengh set to 100bp
    def primerencode(ID, PW):
        IDPW = ID + PW
        res = primerEncoder.str2bin(IDPW)
        seq = primerEncoder.bin2seq(res)
        while len(seq) != 100:
            seq = seq + 'A'
        return seq


    #check No.1 if primer sequence is duplicated in DNA info
    def primercheck_structure(seq,DNA_image,information_DNA,DNA_hos_primer):
        if seq in DNA_image:
            return True
        if seq in information_DNA:
            return True
        if seq in DNA_hos_primer:
            return True
        else:
            return False


    #check No.2 if forward/reverse primer has Tm difference more than 12 degree Celcius
    def primercheck_Tm(seq,DNA_hos_primer):
        count = [0,0,0,0]
        count_2 = [0,0,0,0]
        for i in seq:
            if i == 'A':
                count[0] = count[0]+1
            if i == 'T':
                count[1] = count[1]+1
            if i == 'C':
                count[2] = count[2]+1
            if i == 'G':
                count[3] = count[3]+1
        Tm = 64.9 +41*(count[3]+count[2]-16.4)/(count[0]+count[1]+count[2]+count[3])
        for i in DNA_hos_primer:
            if i == 'A':
                count_2[0] = count_2[0]+1
            if i == 'T':
                count_2[1] = count_2[1]+1
            if i == 'C':
                count_2[2] = count_2[2]+1
            if i == 'G':
                count_2[3] = count_2[3]+1
        Tmp = 64.9 +41*(count_2[3]+count_2[2]-16.4)/(count_2[0]+count_2[1]+count_2[2]+count_2[3])

        if abs(Tm-Tmp) >= 12:
            return False
        else:
            return True


    #check No.3 if primer sequence contains enough GC contents
    def primercheck(seq):
        count = 0
        for i in seq:
            if i == 'G' or i == 'C':
                count = count+1
        GCp = count/len(seq)
        if GCp >= 0.4 and 0.6 >= GCp:
            print("Primer validity successful! You may use this ID/PW combination.")
            return True
        else:
            print("Primer validity failed. please make other combination of ID/PW.")
            return False




class DNAEncoder:
    def __init__(self):
        pass

    def load_dcm(self, in_dir, in_file):
        self.files = []
        self.in_dir = in_dir
        self.in_file = in_file
        
        print("########################## LOADING ##########################")
        self.path = glob.glob(os.path.join(self.in_dir, self.in_file))
        for fname in self.path:
            print("loading: {}".format(fname))
            self.files.append(pydicom.dcmread(fname))
        
        print("file count: {}".format(len(self.files)))
        print("####################### FINISH LOADING ######################")

        # skip files with no SliceLocation (eg scout views)
        self.slices = []
        self.skipcount = 0
        self.non_sorted_slices = []
        for f in self.files:
            if hasattr(f, 'SliceLocation'):
                self.slices.append(f)
            else:
                self.skipcount = self.skipcount + 1
                self.non_sorted_slices.append(f)

        
        print("skipped, no SliceLocation: {}".format(self.skipcount))
        if len(self.slices) == self.skipcount:
            print("No slice location information detected. The image will be stored as the input order.")
            self.length = len(self.files)
            self.slices = self.files
        elif self.skipcount !=0: 
            self.length = len(self.slices) + len(self.non_sorted_slices)
            sort = print("Slice Location information miss. Do you want to sort the files? (yes/no)")
            if sort == 'yes':
                self.slices = sorted(self.slices, key=lambda s: s.SliceLocation)
                self.slices = self.slices + self.non_sorted_slices
            else:
                self.slices = self.slices + self.non_sorted_slices
        elif self.skipcount == 0:
            self.slices = sorted(self.slices, key=lambda s: s.SliceLocation)
            self.length = len(self.slices)


        # ensure they are in the correct order
        self.PixelSpacing = self.slices[0].PixelSpacing
        self.SliceThickness = self.slices[0].SliceThickness
        self.ax_aspect = self.PixelSpacing[1]/self.PixelSpacing[0]
        self.sag_aspect = self.PixelSpacing[1]/self.SliceThickness
        self.cor_aspect = self.SliceThickness/self.PixelSpacing[0]

        # create 3D array
        self.img_shape = list(self.slices[0].pixel_array.shape)
        self.img_shape.append(self.length)
        self.img3d = np.zeros(self.img_shape)
        print(self.img_shape)

        '''
        for i, s in enumerate(self.slices):
            img2d = np.maximum(s.pixel_array,0)
            self.img3d[i, :, :] = img2d
        '''
    def load_hospital_primer(self,primer_file):
        f = open(primer_file, 'r')
        self.DNA_hos_primer = f.read()

    def encode_int(self, x, target_len):
        res = ""
        nuc = ['A', 'T', 'C', 'G']
        for i in range(target_len):
            res = nuc[x & 0b11] + res
            x >>= 2
        return res

    def encode_image(self):
        self.DNA_image = ""
        self.DNA_image += self.encode_int(self.img_shape[0], 8)
        self.DNA_image += self.encode_int(self.img_shape[1], 8)
        self.DNA_image += self.encode_int(self.img_shape[2], 8)

        # fill 3D array with the images from the files
        for i in tqdm(range(self.length)):
            dcm = self.slices[i]
            arr = np.maximum(dcm.pixel_array,0)
            # plt.imshow(arr, cmap=plt.cm.gray)
            # plt.savefig("2.png") 
            img_array = arr.astype('uint16')
            for x in range(self.img_shape[0]):
                for y in range(self.img_shape[1]):
                    self.DNA_image += self.encode_int(img_array[x, y], 8)
            break
        return self.DNA_image

    def encode_patient_info(self):
        information =  str(encoder.files[0]['PatientName'])+'\n'+str(encoder.files[0]['PatientID'])+'\n'+str(encoder.files[0]['PatientBirthDate'])+'\n'+str(encoder.files[0]['PatientAge'])+'\n'+str(encoder.files[0]['PatientSex'])
        ascii_information=[]
        int_list = []
        for i in range(len(information)):
            int_chr = ord(information[i])
            int_list.append(int_chr)
            fourth = np.base_repr(int_chr,base=4,padding=4)
            ascii_information.append(fourth)

        information_DNA = ''
        information_DNA += self.encode_int(len(information),8)
        nuc = ['A', 'T', 'C', 'G']
        for i in range(len(ascii_information)):
            information_i = ascii_information[i]
            for j in range(4):
                information_DNA += nuc[int(information_i[-4+j])]
        
        return information_DNA

    def encode(self,out_file,user_name,password):
        self.user_name = user_name
        self.password = password
        self.DNA = ""
        self.DNA += primerEncoder.primerencode(username, password)
        self.DNA += self.encode_patient_info()
        self.DNA += self.encode_image()
        self.DNA += self.DNA_hos_primer

        f = open(out_file,'w')
        f.write(self.DNA)
        f.close()

class DNADecoder():
    def __init__(self):
        pass

    def decode_int(self, c):
        if c == 'A':
            return 0b00
        elif c == 'T':
            return 0b01
        elif c == 'C':
            return 0b10
        elif c == 'G':
            return 0b11    
    
    def decode_primer(self):
        self.primer = ''

    
    def decode_int_2(self,c):
        if c == 'A':
            return '0'
        elif c == 'T':
            return '1'
        elif c == 'C':
            return '2'
        elif c == 'G':
            return '3'

    def decode_check(self,username, password, dna_file):
        f = open(dna_file, 'r')
        self.DNA_total = f.read()

        self.DNA_primer_input = primerEncoder.primerencode(username, password)
        self.DNA_primer_dc = self.DNA_total[:100]
        if self.DNA_primer_input != self.DNA_primer_dc:
            print("Wrong combination of username and password. Please enter again.")
            return False
        else:
            return True


    def decode_patient_info(self,information_DNA,out_dir,prefix):
        string_information = ''
        for i in range(int(len(information_DNA)/4)):
            DNA_4_i = information_DNA[4*i:4*i+4]
            bin_i = ''
            for j in range(4):
                bin_i += self.decode_int_2(DNA_4_i[j])
            int_i = 0
            for k in range(4):
                int_i+= 4**(3-k)*int(bin_i[k])
            string_information += chr(int_i)
            
        file_name = os.path.join(out_dir, prefix + '.patient_info.txt')
        fo = open(file_name,'w')
        fo.write(string_information)
        fo.close()

    def decode_img(self, DNA_sequence, out_dir, prefix):
        self.prefix = prefix
        DNA_img_sequence = DNA_sequence
        self.x = self.y = self.z = 0
        for i in range(8):
            self.x = (self.x << 2) + self.decode_int(DNA_img_sequence[i])
        for i in range(8):
            self.y = (self.y << 2) + self.decode_int(DNA_img_sequence[i + 8])
        for i in range(8):
            self.z = (self.z << 2) + self.decode_int(DNA_img_sequence[i + 16])
        self.max_id = self.z
        for i in range(self.z):
            file_name = os.path.join(out_dir, self.prefix + '.' + str(i) + '.png')
            img_array = np.zeros((self.x, self.y), dtype='uint16')
            for x in range(self.x):
                for y in range(self.y):
                    img_array[x, y] = 0
                    for k in range(8):
                        img_array[x, y] = (img_array[x, y] << 2) + self.decode_int(DNA_img_sequence[24 + i * self.x * self.y * 8 + (x * self.y + y) * 8 + k])
            plt.imshow(img_array, cmap=plt.cm.gray)
            plt.savefig(file_name) 

    
    def decode(self, dna_file,out_dir,prefix):
        self.prefix = prefix
        f = open(dna_file, 'r')
        self.DNA_total = f.read()

        self.DNA_primer_dc = self.DNA_total[:100]

        self.patient_info_length = 0
        for i in range(100,108):
            self.patient_info_length = (self.patient_info_length << 2) + self.decode_int(self.DNA_total[i]) 
        self.DNA_patient_info_dc = self.DNA_total[100:108+self.patient_info_length*4]
        self.decode_patient_info(self.DNA_patient_info_dc,out_dir,prefix)
        self.DNA_image_dc = self.DNA_total[108+self.patient_info_length*4:-100]

        self.decode_img(self.DNA_image_dc,out_dir,prefix)
        self.DNA_hos_primer_dc = self.DNA_total[-100:]


###################### MAIN ########################
# encoder = DNAEncoder()

# in_dir = input("Please enter the directory: ")
# in_file = input("Please enter the file name (ex.: *.dcm): ")
# primer_file = input("Please enter your hospital primer file (ex.: hospital_primer.txt):")

# # Enter ID and PW
# id_and_pw_check = False
# while id_and_pw_check == False:
#     username = input("ID? (12 characters)")
#     while len(username) >= 13 or 6 >= len(username):
#         print("Your ID should be 7~12 characters.")
#         username = input("ID?")
#     password = input("PW? (12 characters)")
#     while len(password) >= 13 or 6 >= len(password):
#         print("Your PW should be 7~12 characters.")
#         password = input("PW?")

#     seq = primerEncoder.primerencode(username, password)
#     id_and_pw_check = primerEncoder.primercheck(seq)
#     if id_and_pw_check == False:
#         print("Please try another combination of username and password.")


# encoder.load_dcm(in_dir, in_file)
# encoder.load_hospital_primer(primer_file)
# encoder.encode("./dna.txt",username,password)

# # Decode
# decoder = DNADecoder()
# # Enter the DNA file you want to decode
# dna_file = print("Please enter the DNA file you want to decode (ex.: dna.txt):")
# # Enter ID and PW
# id_and_pw_match = False
# while id_and_pw_match == False:
#     username = input("ID? (12 characters)")
#     while len(username) >= 13 or 6 >= len(username):
#         print("Your ID should be 7~12 characters.")
#         username = input("ID?")
#     password = input("PW? (12 characters)")
#     while len(password) >= 13 or 6 >= len(password):
#         print("Your PW should be 7~12 characters.")
#         password = input("PW?")

#     id_and_pw_match = decoder.decode_check(username, password, 'dna.txt')
        

# decoder.decode('dna.txt','./decode_output','result')
# # encoder.encode('./dna.txt')