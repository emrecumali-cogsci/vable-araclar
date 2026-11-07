# Tek Denekli Analiz Aracı — R ile doğrulama / Validation in R

`engine.js`, araçtaki (`../index.html`) hesaplama motorunun aynısıdır. / is the same computation engine embedded in the tool.

1. `Rscript make_ref.R` — SingleCaseES 0.7.4 ve scdhlm (lmeInfo) ile referans değerleri üretir (`ref.json`): 300 rastgele A–B karşılaştırması; scdhlm paketindeki 8 yayımlanmış ve 13 sentetik çoklu başlama düzeyi veri seti (biri boşluklu oturumlu). / produces reference values.
2. `node compare.js` — JS sonuçlarını referanslarla karşılaştırır. / compares JS results with the references.
3. `Rscript ci_tight.R` — NAP güven aralığı köklerini sıkı toleransla (1e-13) yeniden hesaplar. / recomputes NAP CI roots with a tight tolerance.

## Sonuç / Result (2026-11-08)

| Ölçü / Index | En büyük fark / Max. difference |
|---|---|
| NAP, Tau, Tau-U, Tau-BC, PND, PEM, PAND, IRD (tahmin, SE) | < 4·10⁻¹⁵ |
| LRRi, LRRd, SMD (tahmin, SE, GA) | < 6·10⁻¹⁴ |
| NAP/Tau GA (R `uniroot` varsayılan toleransı 1,2·10⁻⁴) | 6·10⁻⁵; sıkı toleransla / with tight tolerance 2,5·10⁻¹¹ |
| BC-SMD g, SE, ν, GA (21 veri seti) | göreli / relative ≤ 1·10⁻⁵ |

Araçtaki "R ile yeniden üret" kodu dört desen için çalıştırıldı ve tablolarla aynı değerleri verdi. / The tool's exported R code was run for all four designs and reproduced the tables.

Formüller SingleCaseES (GPL-3) ve scdhlm/lmeInfo (GPL-3) paketlerindeki yöntemlerle aynıdır; JavaScript'te yeniden yazılmıştır. / Formulas follow the methods in these packages, re-implemented in JavaScript.
