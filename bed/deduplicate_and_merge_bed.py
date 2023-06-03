import sys


def remove_duplicated_lines(lines):
    print('Total lines: ', len(lines))
    lines =  set(lines)
    print('Total not duplicated lines: ', len(lines))
    return lines


def sort_bed_data(bed_data):
    '''The sorted() method sorts tuples by default, using the first item in each tuple.
    '''
    for chr in bed_data.keys():
        print(f'{chr}', len(bed_data[chr]))
        bed_data[chr].sort()


def parse_bed(lines):
    '''Parse bed line to extract position start and end for each chromosome.'''
    bed = {}
    for line in lines:
        chr = line.split('\t')[0]
        start = int(line.split('\t')[1])
        end = int(line.split('\t')[2].strip())
        if chr not in bed:
            bed[chr] = []
        bed[chr].append((start, end))
    sort_bed_data(bed)
    return bed


def merge_intervals(bed_data):
    '''Merge overlapping ranges.
    '''
    bed_data_merged = {}
    for chr, positions in bed_data.items():
        last_start = positions[0][0] # first postion start
        last_end = positions[0][1] # first postion end
        bed_data_merged[chr] = []
        for idx, position in enumerate(positions):
            # current position
            start = position[0]
            end = position[1]
            
            # if current position is equal to the last merged position, skip
            if bed_data_merged[chr]  and (start, end) == bed_data_merged[chr][-1]:
                continue
                                    
            # is the current start in the last position range?
            if start >= last_start and start <= last_end:
                # update the last end with the higher end
                last_end = end if end > last_end else last_end
                # if it is not the last chromossome position
                if  idx != len(positions)  - 1 :
                    continue
            
            # add last position to list because overlaps  are all corrected
            bed_data_merged[chr].append((last_start, last_end))
            
            # update last position with current position and continue
            last_start = start
            last_end = end
        print(f'{chr} total after merge: {len(bed_data_merged[chr])}') 
    return bed_data_merged


def export_merged_bed(bed_file, bed_data):
    '''Export the bed merge data to file
    '''
    print('Exporting the merged bed...')
    with open(bed_file + '.merged', 'w') as f:
        for chr, positions in bed_data.items():
            for position in positions:
                f.write(f'{chr}\t{position[0]}\t{position[1]}\n')


if __name__ == '__main__':
    bed_file = sys.argv[1]
    print(f'Processing {bed_file} ...')
    with open(bed_file, 'r') as f:
        lines = remove_duplicated_lines(f.readlines())
        bed_data = parse_bed(lines)
        bed_data_merged = merge_intervals(bed_data)
        export_merged_bed(bed_file, bed_data_merged)
    print(f'Processing has finished.')
    	
