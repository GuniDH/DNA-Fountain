import random 


# Guni Deyo Haness
# implementation of dna fountain 


class dna_fountain:
    def __init__(self, input, tips_number, tip_deg_map, tip_dna_map):

        # as the input in the instructions comes with space after each group of 4 bits
        self.input=input.replace(' ','')
        # the encoded data, each tip in tips is translated to DNA bases (oligorem)
        self.tips = None
        self.tips_number = tips_number
        # mapping seed values (decimal) to their corresponding tip degrees
        self.tip_deg_map=tip_deg_map
        # mapping binary pair to dna base (encoding) or dna base to binary pair (decoding)
        self.tip_dna_map=tip_dna_map
        # keeping segments as a field so we can restore the random segment choices while decoding
        self.segments=None 


    # input is 32 bit long binary number, we encode it into 16 tips (oligorems)
    def _encode(self):
        # Divide the input into segments of 4 bits
        self.segments = [self.input[i:i+4] for i in range(0, len(self.input), 4)]
        print(f'The 32bit input after dividing to segments in length of 4 bits is {self.segments}')
        # create tips (basically oligomers)
        tips=[]
        for seed in range(self.tips_number):
            random.seed(seed) # setting a seed so the choices can be restored while decoding
            tip_deg=self.tip_deg_map[seed]
            tip_segments = random.sample(self.segments, tip_deg)
            print(f'The segments for the {seed+1}th tip are: {tip_segments}')
            xor_val=0
            for segment in tip_segments:
                xor_val ^= int(segment, 2) # xor is like addition modulo 2
            tip_value = format(seed, '04b') + format(xor_val, '04b') # add seed to beginning of tip value
            # translate each tip to oligomer of dna bases
            oligomer='' # tip after translation to DNA bases is basically an oligomer
            for i in range(0, 8, 2):
                pair = tip_value[i:i+2]
                oligomer += self.tip_dna_map[pair] 
            tips.append(oligomer)
            print(f'The {seed+1}th tip after translation to DNA bases (oligorem) is {oligomer}')
        self.tips=tips


    # decoding oligorems (tips) back to original binary data
    def _decode(self):
        # map tips to their chosen segments
        tip_segments = []
        for oligomer in self.tips:
            # convert to binary from oligomer
            tip_value=''
            for dna_base in oligomer:
                tip_value += self.tip_dna_map[dna_base]
            # each tip value is made out of seed and xor operation between all of its chosen segments (addition modulo 2)
            seed, xor_val = int(tip_value[:4],2), tip_value[4:]
            tip_deg=self.tip_deg_map[seed]
            # regenerate segment choosing (by using same seeds)
            random.seed(seed)
            tip_segment_indexes = random.sample(list(range(8)), tip_deg)
            tip_segments.append([xor_val, tip_segment_indexes])

        # keep track of the indexes of identified segments to know when we finished 
        discovered_segments_indexes=set()
        discovered_segments = len(self.segments)*[" "] 
        while True:
            for _ in range(self.tips_number):
                if len(tip_segments) == 0:
                    break
                # choose random tip
                xor_val, tip_segment_indexes = random.choice(tip_segments)
                if not [xor_val, tip_segment_indexes] in tip_segments:
                    continue # make sure tip wasn't removed
                # we identify a segment once its the only segment of a tip
                if len(tip_segment_indexes)==1:
                    identified_segment_indx = tip_segment_indexes[0]
                    discovered_segments_indexes.add(identified_segment_indx) # keep track of identified segments
                    print(f'Identified segment {identified_segment_indx} with value {xor_val}')
                    # we set the value of the discovered segment to be the xor value of the tip
                    discovered_segments[identified_segment_indx]=xor_val

                    for xv, tsi in tip_segments:
                        if identified_segment_indx in tsi:
                            # removing identified segment from the list of segments of the tip
                            tsi.remove(identified_segment_indx)
                            if tsi:
                                # once we identify a segment whose value is x, we need to go through all its connect tips and perform tip_xor_val = tip_xor_val xor x
                                tip_segments[tip_segments.index([xv, tsi])][0] =  format(int(xor_val,2) ^ int(xv,2), '04b')
                            else: 
                                # in case the identified segment was the only one in the segment list of a tip we totally delete the connection
                                tip_segments.remove([xv, tsi])

            if len(discovered_segments_indexes)==len(self.segments): 
                break

        return ' '.join(discovered_segments)
        

if __name__ == '__main__':
    input='0100 0001 1010 1111 0000 0101 1010 0101'

    tips_number = 16

    tip_deg_map = {
    2: 1, 3: 1, 7: 1, 9: 1, 10: 1, 14: 1,
    0: 2, 1: 2, 4: 2, 6: 2, 11: 2, 13: 2,
    5: 4, 15: 4,
    8: 6,
    12: 7
    }

    tip_dna_map = {
    '00': 'A', '01': 'C', '10': 'G', '11': 'T',
    'A': '00', 'C': '01', 'G': '10', 'T': '11'
    }

    dna_f = dna_fountain(input, tips_number, tip_deg_map, tip_dna_map)
    dna_f._encode()
    output = dna_f._decode()

    print('input:  ', input)
    print('output: ', output)