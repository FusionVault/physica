# Test fixtures

Real ICGEM files, kept whole up to `end_of_head` and truncated to low degrees after it. Downloaded
on 2026-10-23.

| File | Source | Kept | Licence / credit |
|---|---|---|---|
| `costg_2026-07_n10.gfc` | COST-G GRACE-FO July 2026, [ICGEM](https://icgem.gfz.de/getseries/02_COST-G_/Grace-Grace-FO_RL02.1/GSM-2_2026182-2026212_GRFO_COSTG_BF01_0201.gfc) | degree ≤ 10 | CC BY 4.0, doi:10.5880/COST-G.ICGEM_02_L2 |
| `swarm_2025-12_n10.gfc` | Swarm EGF_SHA_2 December 2025, [ESA swarm-diss](https://swarm-diss.eo.esa.int/?do=list&file=swarm%2FLevel2longterm%2FEGF) (`.gfc` inside the ZIP) | degree ≤ 10 | Data provided by the European Space Agency |
| `goco06s_n4.gfc` | GOCO06s, [ICGEM](https://icgem.gfz.de/tom_longtime) | degree ≤ 4 | CC BY 4.0, doi:10.5880/ICGEM.2019.002 |
| `tn14_excerpt.txt` | GRACE TN-14 (GSFC SLR C20/C30, v3, 30 Sep 2026), [PO.DAAC docs](https://archive.podaac.earthdata.nasa.gov/podaac-ops-cumulus-docs/gracefo/open/docs/TN-14_C30_C20_GSFC_SLR.txt) | header, first 3 and last row | NASA, no use constraints |
| `tn13_excerpt.txt` | GRACE TN-13 (CSR RL06.3 geocentre), [PO.DAAC docs](https://archive.podaac.earthdata.nasa.gov/podaac-ops-cumulus-docs/gracefo/open/docs/TN-13_GEOC_CSR_RL0603.txt) | header, first and last month | NASA, no use constraints |
