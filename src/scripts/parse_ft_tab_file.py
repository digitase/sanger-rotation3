#!/usr/bin/env python
#
'''Parse .tab flatfile format
'''

import sys
import errno
import pyparsing as pp
import pandas as pd

class Feature:
    '''A feature record (type=FT) from EMBL flatfile format
    ''' 
    def __init__(self, key=None, loc=None, qualifs=dict()):
        self.type = 'FT'
        self.key = key
        self.loc = loc
        self.qualifs = qualifs
    def from_pp(self, pp_result):
        '''Populate feature from pyparsing output
        '''
        assert(pp_result[0] == self.type), (pp_result[0], self.type)
        self.key = pp_result[1]
        self.loc = pp_result[2]
        if len(pp_result) >= 3:
            self.qualifs = dict((x[0], x[1]) for x in pp_result[3])
    #
    def to_tab(self):
        '''Return string in .tab file format
        '''
        key_col = 6
        loc_col = 22
        type_key_spacer = ' ' * (key_col-len(self.type)-1)
        key_loc_spacer = ' ' * (loc_col-len(self.key)-len(type_key_spacer)-len(self.type)-1)
        type_loc_spacer = ' ' * (loc_col-len(self.type)-1)
        tab_header_template = '{type}{type_key_spacer}{key}{key_loc_spacer}{loc}'
        tab_qualif_template = '{type}{spacer}/{qualif_key}={qualif_value}'
        header = tab_header_template.format(
                type=self.type, 
                type_key_spacer=type_key_spacer,
                key=self.key,
                key_loc_spacer=key_loc_spacer,
                loc=self.loc
        )
        qualifs = [tab_qualif_template.format(
            type=self.type,
            spacer=type_loc_spacer,
            qualif_key=k,
            qualif_value=repr(v) if type(v) is int else '"{}"'.format(repr(v)[1:-1]),
            ) for (k, v) in self.qualifs.items()]
        return '\n'.join([header] + qualifs)
    #
    def __str__(self):
        return(str(map(str, [self.type, self.key, self.loc, self.qualifs])))
    #
    def __repr__(self):
        return(str(map(str, [self.type, self.key, self.loc, self.qualifs])))

class FtTabFile:
    '''A .tab file consisting of multiple FT records
    '''
    @staticmethod
    def get_pp_pattern():
        '''Define parsing pattern for FT record type
        '''
        ft_type = pp.Literal('FT')
        ft_key = pp.Word(pp.printables)
        ft_loc = pp.Word(pp.nums + '.')
        #
        ft_header = ft_type + ft_key + ft_loc
        #
        ft_qualif_key = pp.Word(pp.printables.replace('=' ,''))
        ft_qualif_value = pp.restOfLine()
        ft_qualif = pp.Suppress(ft_type) + pp.Group(pp.Suppress('/') + ft_qualif_key + pp.Suppress('=') + ft_qualif_value)
        #
        ft = pp.Group(ft_header + pp.Group(pp.ZeroOrMore(ft_qualif)))
        pattern = pp.ZeroOrMore(ft)
        return(pattern)
    #
    def __init__(self):
        self.type = "FT"
        self.features = []
    #
    def parse_file(self, filename):
        '''Apply pyparsing pattern to .tabfile
        '''
        pattern = FtTabFile.get_pp_pattern()
        result = pattern.parseFile(filename)
        for x in result:
            ft = Feature()
            ft.from_pp(x)
            self.features.append(ft)
        print("Parsed {} records.".format(len(self.features)))
    def __str__(self):
        return(str(self.features))
    #
    def write_csv(self, out_file=None):
        '''Write features to .csv
        '''
        out_df = pd.DataFrame()
        for k in ['type', 'key', 'loc']:
            col = [getattr(f, k) for f in self.features]
            out_df.loc[:, k] = col
        # Get set of keys in feature qualifiers
        all_qualif_keys = set()
        for f in self.features:
            all_qualif_keys.update(f.qualifs.keys())
        all_qualif_keys = list(all_qualif_keys)
        for k in all_qualif_keys:
            col = [f.qualifs[k] if k in f.qualifs else None for f in self.features]
            out_df.loc[:, k] = col
        #  Write out csv
        return(out_df.to_csv(out_file))
    #
    def write_tab(self, out_file=None):
        '''Write features to .tab
        '''
        output = '\n'.join(f.to_tab() for f in self.features)
        if out_file:
            with open(out_file, 'w') as out_fhandle:
                out_fhandle.write(output)
        else:
            print(output)
        
if __name__ == "__main__":
    test_ft_tab_string = '''FT   SNP             2654089
FT                   /strand="reverse"
FT                   /node="2->7521_5#31"
FT                   /SNP="T->G"
FT                   /codon_type="Nonsynonymous"
FT                   /colour=4
FT                   /taxa="7521_5#31"
FT   SNP             1936327
FT                   /node="836->838"
FT                   /SNP="T->G"
FT                   /homoplasy="convergence with branch leading to 7469_7#82, convergence with branch leading to 242, convergence with branch leading to 7480_7#17, convergence with branch leading to 7480_8#27, convergence with branch leading to 7554_6#82, convergence with branch leading to 12641_1#70, convergence with branch leading to 7480_8#41, convergence with branch leading to 7414_8#64, convergence with branch leading to 7480_8#10, convergence with branch leading to 7480_8#2, convergence with branch leading to 7554_6#34, convergence with branch leading to 7554_6#18"
FT                   /codon_type="Intergenic"
FT                   /colour=2
FT                   /taxa="7712_8#5, 7480_7#23"'''
    
    in_file = sys.argv[1]
    out_file = sys.argv[2]

    try:
        foo = FtTabFile()
        foo.parse_file(in_file)
        foo.write_csv(out_file)
    except IOError as e:
        print("({})".format(e))

