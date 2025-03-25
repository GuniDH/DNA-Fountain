import random 
import numpy as np
from collections import defaultdict

# Author: Guni 
# Improved DNA fountain

class dna_fountain:
    def __init__(self, tip_deg_map, tip_dna_map, letters_per_word, bits_per_word, bits_for_barcode, bits_for_seed, bits_for_deg, bits_for_segment, tip_copies_num, tips_num, segment_num,orig_input_size):
        self.tip_deg_map=tip_deg_map
        self.tip_dna_map=tip_dna_map
        self.letters_per_word=letters_per_word 
        self.bits_per_word=bits_per_word 
        self.bits_for_barcode=bits_for_barcode
        self.bits_for_seed=bits_for_seed
        self.bits_for_deg=bits_for_deg
        self.bits_for_segment = bits_for_segment
        self.tip_copies_num=tip_copies_num
        self.tips_num = tips_num
        self.segment_num = segment_num
        self.orig_input_size=orig_input_size

    def _get_barcode(self, barcodes): # generate unique barcode
        barcode = None
        while True:
            barcode=random.randint(0, (2**self.bits_for_barcode) - 1)  # random bits_for_barcode-bit barcode
            if barcode not in barcodes:
                barcodes.add(barcode)
                break  
        return barcode 

    def _encode(self, input):
        # Divide the input into segments 
        segments = [input[i:i+self.bits_for_segment] for i in range(0, len(input), self.bits_for_segment)]
        print(f'The input after dividing to segments is {segments}')

        oligomers=[]
        barcodes = set()
        for seed in range(self.tips_num):
            random.seed(seed) # setting a seed so the choices can be restored while decoding
            tip_deg=self.tip_deg_map[seed]
            tip_segments = random.sample(segments, tip_deg)
            print(f'The segments for the {seed+1}th tip are: {tip_segments}')
            xor_val=0
            for segment in tip_segments:
                xor_val ^= int(segment, 2) # xor is like addition modulo 2
            barcode = self._get_barcode(barcodes)
            payload=xor_val
            tip_value = format(barcode, f'0{self.bits_for_barcode}b') + format(seed, f'0{self.bits_for_seed}b') + format(tip_deg, f'0{self.bits_for_deg}b') + format(payload, f'0{self.bits_for_segment}b') 
            # translate each tip to oligomer of dna bases
            oligomer='' 
            for i in range(0, self.bits_for_barcode+self.bits_for_seed+self.bits_for_deg+self.bits_for_segment, self.bits_per_word): 
                trio = tip_value[i:i+self.bits_per_word]
                word = self.tip_dna_map[trio]
                word = word[0] + self.tip_dna_map[word[1]] # replace M or K by 2 random bases from [A,C,G,T]
                oligomer += word
            print(f'The {seed+1}th tip after translation to DNA bases (oligomer) is {oligomer}')
            for _ in range(self.tip_copies_num): # create self.tip_copies_num (100) sequences for each tip
                oligomers.append(oligomer)
        random.shuffle(oligomers) # This step represents the mixing of all the DNA sequences in one test tube containing the entire coded file for storage
        return oligomers

    # decoding oligomers (tips) back to original binary data
    def _decode(self, oligomers):
        # map tips to their chosen segments
        tip_segments = []
        barcode_to_tips_map = defaultdict(list)
        for oligomer in oligomers:
            # convert to binary from oligomer
            tip_value=''
            for i in range(0, len(oligomer), self.bits_per_word):
                first_letter, second_letter = oligomer[i], self.tip_dna_map[oligomer[i+1:i+1+self.letters_per_word]]
                word = first_letter+second_letter
                tip_value += self.tip_dna_map[word]
            # each tip value is made out of barcode, seed, degree and payload - xor operation between all of its chosen segments (addition modulo 2)
            barcode = int(tip_value[:self.bits_for_barcode],2)
            seed = int(tip_value[self.bits_for_barcode:self.bits_for_barcode+self.bits_for_seed],2)
            tip_deg = int(tip_value[self.bits_for_barcode+self.bits_for_seed:self.bits_for_barcode+self.bits_for_seed+self.bits_for_deg],2)
            payload = tip_value[self.bits_for_barcode+self.bits_for_seed+self.bits_for_deg:]
            # regenerate segment choosing (by using same seeds)
            random.seed(seed)
            tip_segment_indexes = random.sample(list(range(self.segment_num)), tip_deg)
            tip_segments.append([barcode,payload, tip_segment_indexes])

        # keep track of the indexes of identified segments to know when we finished 
        discovered_segments_indexes=set()
        discovered_segments = self.segment_num*[" "] 
        while True:
            for _ in range(self.tips_num):
                if len(tip_segments) == 0:
                    break
                # choose random tip
                barcode, payload, tip_segment_indexes = random.choice(tip_segments)
                if not [barcode, payload, tip_segment_indexes] in tip_segments:
                    continue # make sure tip wasn't removed
                # we identify a segment once its the only segment of a tip
                identified_segment_indx = tip_segment_indexes[0]
                if len(tip_segment_indexes)==1 and identified_segment_indx not in discovered_segments_indexes:
                    barcode_to_tips_map[barcode].append(payload)
                    discovered_segments_indexes.add(identified_segment_indx) # keep track of identified segments
                    print(f'Identified segment {identified_segment_indx} with value {payload}')
                    # we set the value of the discovered segment to be the xor value of the tip (payload)
                    discovered_segments[identified_segment_indx]=payload
                    for b, pl, tsi in tip_segments:
                        if identified_segment_indx in tsi:
                            # removing identified segment from the list of segments of the tip
                            tsi.remove(identified_segment_indx)
                            if tsi:
                                # once we identify a segment whose value is x, we need to go through all its connect tips and perform tip_xor_val = tip_xor_val xor x
                                tip_segments[tip_segments.index([b, pl, tsi])][1] =  format(int(payload,2) ^ int(pl,2), f'0{self.bits_for_segment}b')
                            else: 
                                # in case the identified segment was the only one in the segment list of a tip we totally delete the connection
                                tip_segments.remove([b, pl, tsi])

            if len(discovered_segments_indexes)==self.segment_num: # found all segments
                break
        return (''.join(discovered_segments))[:self.orig_input_size], barcode_to_tips_map # remove padding
    
    @staticmethod
    def ideal_soliton_distribution(k):
        """Generate the ideal soliton distribution."""
        rho = np.zeros(k + 1)
        rho[1] = 1 / k
        for d in range(2, k + 1):
            rho[d] = 1 / (d * (d - 1))
        return rho

    @staticmethod
    def robust_soliton_distribution(k, c=0.1, delta=0.5):
        """Generate the robust soliton distribution."""
        R = c * np.log(k / delta) * np.sqrt(k)
        tau = np.zeros(k + 1)
        for d in range(1, k + 1):
            if d < k / R:
                tau[d] = R / (d * k)
            elif d == k / R:
                tau[d] = R * np.log(R / delta) / k
            else:
                tau[d] = 0
        
        rho = dna_fountain.ideal_soliton_distribution(k)
        mu = rho + tau
        mu /= mu.sum()  # Normalize to make it a probability distribution
        return mu

    @staticmethod
    def generate_tip_degrees(num_tips, max_degree):
        """Generate a dictionary mapping tip indices to degrees using the robust soliton distribution and scale to max_degree with minimum degree 1."""
        mu = dna_fountain.robust_soliton_distribution(num_tips)
        degrees = np.random.choice(np.arange(1, num_tips + 1), size=num_tips, p=mu[1:])
        # Scale the degrees so that the minimum is 1 and the maximum is max_degree
        min_generated_degree = min(degrees)
        max_generated_degree = max(degrees)
        # Linear scaling to ensure degrees range from 1 to max_degree
        scaled_degrees = 1 + ((degrees - min_generated_degree) * (max_degree - 1)) / (max_generated_degree - min_generated_degree)
        # Convert to integer values
        scaled_degrees = np.floor(scaled_degrees).astype(int)
        # Ensure that the scaled degrees do not exceed max_degree
        scaled_degrees = np.minimum(scaled_degrees, max_degree)
        # Create the tip degree mapping
        tip_deg_map = {}
        for i, deg in enumerate(sorted(scaled_degrees)):
            tip_deg_map[i] = int(deg)
        return tip_deg_map

if __name__ == '__main__':
    
    input='000111111000101100001000000010000010111000010010100111100110000100000000000000110011000000110000001100010010111001101010011100000110011100000000001011000111101001111011001110001101001111111111111111111111111011100100100110000100001001101010110010100110000111001000001000010011011010001010001100000110010111001110111000110110110110011000111000111110011010011000110100111011101100010000010001100110010011100110001011000111001000101000011010101100110100010000010100101100101100101001110110110001110001100010100100111001000001000011100101001101001000011100111010100011110110000111001110010001010100010001101000110000100100100001001110011110011011010100011011111001111111101111111101010111101110111101111111100111110000111101101011111110011111110101011110000011111010011110111101110111110100111111111011101111101110111010010111100111111100111111111111111001110100000110100111000100000000100010101011000001000100000000000111100001111000000000100000001000011111111011000000101111111001001110000000001100110000000001110001111000010010000100100001001000010100000100100011110000100100011111000101010001011000111001011111100100110001000100000111000111010001010010010111000100110001001100010111000101000101001010000110100000100010000010001010001010101110101011100000011001010101010101010101010010111101101000100110110101111010111110011100001101000111101000101000101010101010101010100111100010110111001100110010001101110011000010110010101101101001001010111111011001001010000011100010110011110111000010110110011101010011010010110010101110001001101111000001110111011100001011100100010111010011111110110011110011110000111100111000101000000010111111000111110000000011110010110111001100110111111001011110011111111011101111000000011111110011111111100001111110011100001000001011110001111110111111100000001010000011010000010100011110000101100011110111000000010110001101001001110000000001000111011110000111100001111001100011111111000011110001'
    if not input or input == '':
        raise Exception('Input cannot be empty')
    bits_for_segment = min(len(input),24) # 24 divisible by bits_per_word = 3
    padded_input = input + '0'*(bits_for_segment-(len(input)%bits_for_segment)) if len(input)%bits_for_segment !=0 else input # pad input with zeros to make it divisible by segment size
    segment_num = len(padded_input)//bits_for_segment # amount of segments
    tips_num = 4*segment_num # in mamman 15 the tips amount was twice segment_num, but the input there was very small.
                                # after performing many tests on many inputs with different sizes (even bigger than the input here) I found 4 times the size to always work.
    tip_deg_map = dna_fountain.generate_tip_degrees(tips_num,segment_num)  # mapping seed values to their corresponding tip degrees
    bases=['A','C','G','T'] # extra bases are M and K such that each will be represented by 2 random bases
    m_bases=random.sample(bases, 2)
    k_bases=list(set(bases) - set(m_bases)) # no need to choose 2 random bases because the whole set is just 4 items so only two are left after the first choice
    m_map, k_map=''.join(m_bases),''.join(k_bases) # map M,K to two different random bases, each
    tip_dna_map = { # language representation 
    '000': 'AM', '010': 'CM', '100': 'GM', '110': 'TM',
    'AM': '000', 'CM': '010', 'GM': '100', 'TM': '110',
    '001': 'AK', '011': 'CK', '101': 'GK', '111': 'TK',
    'AK': '001', 'CK': '011', 'GK': '101', 'TK': '111',
    'M': m_map, 'K': k_map, m_map: 'M', k_map: 'K'
    }
    letters_per_word=2
    bits_per_word=3 
    bits_for_barcode = 10
    # -2 to ignore 0b.
    bits_for_deg = len(bin(segment_num))-2  # the max value for a degree is the amount of segments
    bits_for_seed = len(bin(tips_num-1))-2 # theres a seed for each tip and seed's value is the index (starting from 0) so the maximum val is tips_num-1
    s = bits_for_barcode+bits_for_seed+bits_for_deg+bits_for_segment
    bits_for_seed = bits_for_seed + (60-(s%60)) if s%60 !=0 else bits_for_seed # I chose packet size to be 60 (divisible by bits_per_word and most efficient after performing grid search), so i make sure each packet is divisible by 60
    tip_copies_num=100

    dna_f = dna_fountain(tip_deg_map, tip_dna_map, letters_per_word, bits_per_word, bits_for_barcode, bits_for_seed, bits_for_deg, bits_for_segment, tip_copies_num, tips_num, segment_num,len(input))
    oligomers = dna_f._encode(padded_input)
    output, barcode_to_tips_map = dna_f._decode(oligomers)
    print('barcode to tips mapping: ', barcode_to_tips_map)
    print('input:  ', input)
    print('output: ', output)
    print(f'input is {'' if input==output else 'not '}equal to output')
