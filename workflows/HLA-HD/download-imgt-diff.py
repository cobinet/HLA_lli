import pandas as pd
import requests


class VersionGenerator:
    def __init__(self, start=32, end=50) -> None:
        self.versions = range(start, end + 1)

    def __iter__(self):
        self.iterator = iter(self.versions)
        return self

    def __next__(self):
        return f"v3{next(self.iterator)}"


total_mod = pd.DataFrame()
total_add = pd.DataFrame()
for v in VersionGenerator():
    print("downloading version:", v)
    url = f"https://www.ebi.ac.uk/ipd/imgt/hla/release/{v}/"
    html = requests.get(url).content
    try:
        df_list = pd.read_html(html)
        added = df_list[0]
        modified = df_list[2]
        modified["ver"] = v
        added["ver"] = v
        total_mod = pd.concat([total_mod, modified])
        total_add = pd.concat([total_add, added])
    except:
        continue
total_mod.to_csv("modified.csv")
total_add.to_csv("added.csv")
