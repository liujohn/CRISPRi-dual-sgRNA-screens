def getMismatchDict(table, column, trim_range=None, allowOneMismatch=True):
    mismatch_dict = dict()

    if trim_range:
        col = table[column].apply(lambda seq: seq[trim_range[0]:trim_range[1]])
    else:
        col = table[column]

    for sgRNA, seq in col.iteritems():
        if seq in mismatch_dict:
            print('clash with 0 mismatches', sgRNA, seq)
            mismatch_dict[seq] = 'multiple'
        else:
            mismatch_dict[seq] = sgRNA

        if allowOneMismatch:
            for position in range(len(seq)):
                mismatchSeq = seq[:position] + 'N' + seq[position + 1:]

                if mismatchSeq in mismatch_dict:
                    print('clash with 1 mismatch', sgRNA, mismatchSeq)
                    mismatch_dict[seq] = 'multiple'
                else:
                    mismatch_dict[mismatchSeq] = sgRNA

    return mismatch_dict


def matchBarcode(mismatch_dict, barcode, allowOneMismatch=True):
    if barcode in mismatch_dict:
        match = mismatch_dict[barcode]

    elif allowOneMismatch:
        match = 'none'
        for position in range(len(barcode)):
            mismatchSeq = barcode[:position] + 'N' + barcode[position + 1:]

            if mismatchSeq in mismatch_dict:
                if match == 'none':
                    match = mismatch_dict[mismatchSeq]
                else:
                    match = 'multiple'

    else:
        match = 'none'

    return match
