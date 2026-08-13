# 2026 Bundibugyo Ebola Outbreak (DRC/Uganda) — Data for tsgc analysis

Source: Institut National de Sante Publique (INSP), DRC — SitRep MVE series,
transcribed and harmonised by INRB-UMIE/BDBV2026-Data
(https://github.com/INRB-UMIE/BDBV2026-Data), retrieved 2026-08-13.

Reuse: attribute to INSP and cite the specific sitrep number/date; confirm
distribution terms with INSP (pierre.akilimali@insp.cd) before external
republication. Cite INRB-UMIE/BDBV2026-Data if reusing their processing.

Context: WHO declared this outbreak (Bundibugyo ebolavirus) a Public Health
Emergency of International Concern on 16 May 2026. As of the most recent
national file here (2026-08-10), DRC had reported 4,449 confirmed cases and
2,061 deaths, making it the second-largest Ebola outbreak ever recorded
(after West Africa 2014-16) and the fastest-growing on record.

## Files

National-level (nom = "DRC"), one row per report date:
- insp_sitrep__national_cumulative_confirmed_cases__daily.csv   (2026-05-14 to 2026-08-10)
- insp_sitrep__national_cumulative_confirmed_deaths__daily.csv  (2026-05-14 to 2026-08-10)
- insp_sitrep__national_cumulative_suspected_cases__daily.csv
- insp_sitrep__national_cumulative_suspected_deaths__daily.csv
- insp_sitrep__national_cumulative_recovered_cases__daily.csv

Health-zone level (nom = health zone name, e.g. "Mongbwalu", "Bunia", "Rwampara"),
61 zones, dates from 2026-05-14 to 2026-08-04:
- zone_cumulative_confirmed_cases__daily.csv
- zone_cumulative_confirmed_deaths__daily.csv
- zone_new_confirmed_cases__daily.csv   (daily incidence by zone, sparsely reported)

Note: reporting is irregular (not every calendar day has a sitrep), and
occasional downward revisions occur when a later report corrects an earlier
cumulative total. Both are relevant when setting up idx_series/idx_calendar
in tsgc (irregular reporting dates -> non-daily/patterned calendar; revisions
-> decide whether to model as reported-at-the-time or use latest-vintage
values).

SOURCE_metadata.yaml is the original metadata.yaml from the source repo
(citation, license, retrieval date, processing notes).
