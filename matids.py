import dotenv
import pandas as pd
from mp_api.client import MPRester

m = MPRester(dotenv.get_key(".env", "MATERIAL_PROJ_API_KEY"))
data_one = m.materials.get_data_by_id("mp-1030")
print(data_one)
matid = data_one.material_id
nsites = data_one.nsites
density = data_one.density
volume = data_one.volume
elements = data_one.elements
comp = data_one.composition
sym = data_one.symmetry

df = pd.DataFrame ({
    "material_id": [matid],
    "nsites": [nsites],
    "density": [density],
    "volume": [volume],
    "elements": [elements],
    "composition": [comp],
    "symmetry": [sym]
})
df.to_excel("mp-1030.xlsx", index=False)