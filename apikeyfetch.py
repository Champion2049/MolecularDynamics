import dotenv
from mp_api.client import MPRester

m = MPRester(dotenv.get_key(".env", "MATERIAL_PROJ_API_KEY"))
data_one = m.materials.get_data_by_id("mp-1030")
print(data_one)