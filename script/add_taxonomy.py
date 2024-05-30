import sys
import pandas as pd

tax = sys.argv[1]
abu = sys.argv[2]
output_name = sys.argv[3]

taxinfo = pd.read_csv(tax, sep="\t", header=0)
abundance = pd.read_csv(abu, sep="\t", header=0)

taxinfo = taxinfo[["VCID", "VC_linegae","VC_lifestyle","VC_host_lineage"]]
#taxinfo.set_index('VCID', inplace=True)
#taxinfo_dict = taxinfo.to_dict('index')
#taxinfo_dict = {index: row.tolist() for index, row in taxinfo.iterrows()}

abundance = abundance[["VCID", "Relative_abudannce"]]


taxinfo['VCID'] = taxinfo['VCID'].astype(str)
abundance['VCID'] = abundance['VCID'].astype(str)


abundance = pd.merge(abundance, taxinfo, on='VCID', how='left')
abundance = abundance.rename(columns={'Relative_abudannce':'Relative_abundance', 'VC_linegae':'VC_lineage'})


abundance.to_csv(output_name,sep="\t", index=False)