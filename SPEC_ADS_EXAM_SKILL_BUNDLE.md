---
title: ADS Online Problem-Solving Exam Skill Bundle Spec
status: draft
owner: Eric Yiru Chen
created: 2026-05-20
implementation_status: not_started
review_status: second_pass_completed
---

# ADS Online Problem-Solving Exam Skill Bundle Spec

## Objective

Create a Codex/agent skill bundle that helps solve ADS/ADS2 online exam questions using the user's local course materials as the source of truth. The bundle should make an AI agent faster and more reliable at recognizing question types, selecting the correct statistical method, writing correct R code, interpreting output, and formatting final exam answers.

The bundle is not intended to copy entire course documents into prompts. It should distill reusable workflows, decision rules, code patterns, interpretation templates, and common exam traps from the provided materials.

## Source Material

Primary folders:

- `/Users/eric_yiru/Library/CloudStorage/OneDrive-InternationalCampus,ZhejiangUniversity/过去的/大二下/incourse/ADS2`
- `/Users/eric_yiru/Library/CloudStorage/OneDrive-InternationalCampus,ZhejiangUniversity/过去的/大二下/incourse/大二下`

Important observed file groups:

- ADS2 review and mock material:
  - `ADS2_CC1_2023-24_Review.pdf`
  - `ADS2_mock_CC2_2023-24.pdf`
  - `mock.Rmd`
  - `mock.pdf`
  - `ADS2Week2.4_Practical.pdf`
  - `Categorical_data(2).pdf`
  - `2023_ADS2_S2W10L1_preview.pdf`
  - `ADS2_W24L1_preview.pdf`
- Past ADS/ADS2 papers and scripts:
  - `19级ADS2 S2.pdf`
  - `20级ADS2 S2.pdf`
  - `19ADS.Rmd`
  - `20ADS.Rmd`
- Practical RMarkdown files:
  - `ADS2_CorLR_Practical 2023-24(1).Rmd`
  - `Correlation and regression.R`
  - `regression.Rmd`
  - `week4 pra.Rmd`
  - `week7.Rmd`
  - `week10.Rmd`
  - `week14.Rmd`
  - `feature_extract.Rmd`
- Lecture slides:
  - `week2.4_pre_lecture_slides.pptx`
  - `week2.4_lecture_slides.pptx`
  - `week2.13_lecture_slides(1).pptx`
  - `Lecture2.14_slides.pptx`
- Datasets likely used in examples or exams:
  - `WT.csv`
  - `WBpopulation.csv`
  - `KO.csv`
  - `mnist_train.csv`
  - `really_tiny_dataset.csv`
  - `owid-covid-data.txt`
  - `t1d_drug.csv`
  - `guests.csv`
  - `ADS_files/*.csv`
  - `ICA/*.csv`

Non-ADS files in the second folder, such as IFBS, GP2, and DST2 PDFs, should be ignored unless they are explicitly referenced by ADS files.

Second-pass review sources added on 2026-05-20:

- User-provided exam notes in this conversation, especially:
  - data cleaning template: `str`, `head`, `colSums(is.na())`, `na.omit`, `duplicated`, `unique`, typo inspection, `gsub`, long/wide reshaping;
  - Bayes/turtles worked example;
  - k-means neuron classification workflow;
  - power-analysis writeup templates for chi-square, ANOVA, and t-tests;
  - categorical-data decision notes;
  - t-test, Wilcoxon, ANOVA, and report-writing templates.

Extraction notes from the second review pass:

- PDF text extraction succeeded for the main ADS2 PDFs using `pdftotext`.
- PPTX text extraction succeeded for `Lecture2.14_slides.pptx`, `week2.4_lecture_slides.pptx`, and `week2.4_pre_lecture_slides.pptx`.
- `week2.13_lecture_slides(1).pptx` is 0 bytes and should be marked unusable in the source index.
- The relevant past/mock exam PDFs have consistent structure: three questions, RMarkdown/PDF submission, code plus explanatory text, method choice, hypotheses, assumptions, results, discussion, and "what next" recommendation.

## Review Findings: Actual Common Exam Tasks

The skill bundle should prioritize the following task families because they recur in past papers, mock papers, ICA work, or course practicals.

### 1. Exam Report Workflow

Every coding challenge expects a knitted RMarkdown-style response, not just a numeric answer. The skill must force this answer shape:

1. import and inspect data;
2. clean or reshape only as needed;
3. make a useful plot;
4. state the method and why it fits;
5. state null and alternative hypotheses;
6. check assumptions or justify non-parametric/simulation alternative;
7. run the test/model;
8. interpret the p-value/effect/estimate in context;
9. make a practical recommendation and mention limitations.

This workflow is as important as choosing the right test.

### 2. Data Cleaning and Reshaping

Common exam expectation:

- Load with `read.csv()` or `read.table()` from the working directory.
- Inspect with `str()`, `head()`, `summary()`, `colSums(is.na())`, `anyNA()`, `duplicated()`, and `unique()`.
- Remove missing rows only if justified: `na.omit()` or targeted filtering.
- Remove exact duplicates with `unique()` or `!duplicated()`, but only after checking whether duplicates are true data-entry duplicates.
- Inspect text/factor columns with `unique()`; fix obvious typos with targeted `gsub()` or recoding.
- Convert categorical predictors to `factor`, especially dose/group variables before ANOVA or plotting.
- Reshape when the question needs paired before/after data, repeated measures, or tidy plotting.

Required recipes:

- `pivot_longer()` / `pivot_wider()` as the modern default.
- Include `gather()` / `spread()` only as legacy course-compatible alternatives.
- Combine units when needed, for example minutes plus seconds into one time variable.
- Merge before/after files by patient/sample ID, then check pairing.

### 3. Probability, Bayes, and Conditional Reasoning

Observed tasks:

- Turtles/beaches problem: compute posterior probability using Bayes' theorem and state prior assumptions.
- Lie-detector problem: conditional probability with imperfect sensitivity/specificity and base-rate effects.
- Card-problem logic: identify which observations can falsify an implication.
- Markov-chain simulation: draw states/transitions and simulate path length or hitting time.
- Bayes factor: compare two hypotheses using likelihood ratio times prior ratio.

Skill requirements:

- Make priors explicit.
- Define event notation before computing.
- Distinguish `P(A|B)` from `P(B|A)`.
- Use probability trees or formulas when the problem is simple.
- Use simulation when `P(D|H)` or hitting time is awkward to derive analytically.

### 4. Power and Sample Size

Observed tasks:

- One-sample t-test power by simulation and `power.t.test`.
- Two-sample and paired sample-size comparison for diet-pill style problems.
- Chi-square power follow-up when small cells or future sample-size recommendation appears.
- Course emphasis: recommendations to increase sample size should only get credit when supported by power analysis.

Mathematical points:

- Power is `1 - beta`, the probability of rejecting a false null hypothesis.
- Low power means high Type II error risk; non-significance does not prove no effect.
- Current power and required sample size should be reported with alpha, effect size, and test type.
- For `pwr.t.test`, `n` is per-group sample size for two-sample tests, not total sample size.
- For paired designs, use `type = "paired"` or simulation of paired differences.
- For ANOVA power, convert eta-squared to Cohen's `f` as `sqrt(eta_sq / (1 - eta_sq))`.
- For chi-square power, use a valid effect size `w`; for contingency tables this is commonly Cramer's-V-like: `sqrt(X2 / (N * min(r - 1, c - 1)))`, while goodness-of-fit can use `sqrt(X2 / N)` or the probability-vector formula.

### 5. t-Tests and Wilcoxon Tests

Observed tasks:

- Paired before/after running-time data from CC1 review.
- Independent Monday/Sunday or male/female comparisons where normality may fail.
- One-sample t-test and simulation tasks in the power practical.
- Wilcoxon signed-rank test for paired ordinal/non-normal SCI before/after data.

Decision rules:

- Paired data: test normality of the differences, not each raw group separately.
- Independent two-group numeric data: check normality and independence; use Welch t-test by default if variances are uncertain.
- If non-normal or ordinal:
  - independent groups: Wilcoxon rank-sum / Mann-Whitney U, `wilcox.test(x, y, paired = FALSE)`;
  - paired groups: Wilcoxon signed-rank, `wilcox.test(x, y, paired = TRUE)`.
- One-tailed alternatives must match the question wording and the order of arguments.

### 6. ANOVA and Factorial Designs

Observed tasks:

- Mock vitamin C/tooth-growth question with supplement and dose factors.
- Course review notes on assumptions and post-hoc tests.
- User notes currently say "only one-way"; this must be corrected. The mock exam uses a two-factor design, and the skill should support one-way and two-way/factorial ANOVA with interaction.

Decision rules:

- Use ANOVA when the response is numeric and predictors are categorical groups/factors.
- Use `aov(y ~ group)` for one factor.
- Use `aov(y ~ factor1 * factor2)` when there are two factors and interaction may matter.
- Check:
  - independent sampling from the study design;
  - residual normality: `shapiro.test(resid(model))` and/or residual histogram/QQ plot;
  - approximate equality of variance: residuals vs fitted plot, optionally Levene test if available.
- Use `TukeyHSD(model)` after significant ANOVA to identify which groups differ.
- Interpret main effects carefully if interaction is significant.

### 7. Categorical Data and Chi-Square/Fisher Tests

Observed tasks:

- Goodness-of-fit: season preferences, Mendelian genotype ratios, expected 1:2:1 or sex-genotype combinations.
- Independence/homogeneity: allergy vs season, gene knockout survival, satisfaction by opening time, dementia/footballer count.
- Small expected cells: check chi-square assumptions and use Fisher's exact test when inappropriate.
- 3-way categorical examples: gene, sex, survival.

Decision rules:

- Goodness-of-fit: one categorical variable vs expected probabilities; use `chisq.test(x, p = expected_probs)`.
- Independence/homogeneity: contingency table of two categorical variables; use `chisq.test(table(x, y))`.
- Check expected counts: no expected count below 1 and no more than 20% below 5.
- Use `fisher.test()` when expected counts are too small, especially for small contingency tables.
- For ordinal categorical variables, do not treat ranks as interval-scale numeric without justification; consider Wilcoxon or Kruskal-Wallis.
- `CrossTable()` is useful for display, but run `chisq.test()` on the actual table/matrix, not on the printed object.

### 8. Correlation and Regression

Observed tasks:

- COVID vaccination vs new cases: correlation, simple regression, and interpretation without causal overclaiming.
- Local/time-window regression to detect increase/decrease periods.
- Multiple regression with additional predictors and interactions.
- Comparing slopes/correlations across groups.
- ICA opioid trend: fit per-age-group regressions, compare slopes between time periods.

Decision rules:

- Correlation answers association, not causation.
- Regression requires numeric response and appropriately encoded predictors.
- Linear regression assumptions:
  - independent errors;
  - linear relationship;
  - homoscedastic residuals;
  - approximately normal residuals;
  - no severe multicollinearity for multiple regression.
- Use `summary(lm_fit)` for coefficients, p-values, R-squared.
- Use `plot(lm_fit)` for residual diagnostics.
- Use `anova(model_simple, model_complex)` to compare nested regression models.
- Use VIF if available for multicollinearity; otherwise inspect predictor correlations and model instability.

### 9. Clustering and Feature Extraction

Observed tasks:

- Neuron classification with `vmndata.csv`: plot original labels, apply `kmeans`, compare clusters to original classification, test subsets of fit parameters.
- Feature extraction from tiny image/MNIST-like data: reshape pixels, compute row/column summaries, train simple classifier in practice files.

Skill scope:

- Include k-means clustering as an exam task because it appears in a past ADS2 paper.
- Include feature extraction as a lower-priority course-support recipe, not a first-priority exam method unless a question explicitly asks classification/features.
- For k-means:
  - set seed;
  - scale features when units differ;
  - use `nstart`;
  - plot clusters;
  - compare to known labels with a table and, if available, adjusted Rand index.
- `factoextra` and `CommKern` may not be installed; provide base-R fallback code.

### 10. Bootstrapping and Simulation

Observed tasks:

- Bootstrapping lecture: use sample-with-replacement when assumptions fail or distribution is unknown.
- Week14 notes: bootstrap medians and proportions.
- Power practical: simulate repeated experiments.
- Coffee shop answer used bootstrap confidence intervals for satisfaction proportions.

Skill requirements:

- Provide a generic simulation skeleton:
  - set seed;
  - define statistic;
  - repeat many times;
  - compute p-value or confidence interval;
  - interpret in context.
- Distinguish bootstrap confidence intervals from null-randomization tests.

### 11. ICA/Open Analysis

Observed tasks:

- `substance_use.csv` ICA:
  - filter by measure/location/sex/age/cause/year;
  - summarize `val`, `upper`, `lower`;
  - visualize trends over time and by age/sex;
  - answer specific questions;
  - ask a defensible original question;
  - interpret clearly and avoid overambitious unsupported modeling.
- ICA grading emphasizes clarity, complete code, interpretation, reproducibility, and code that would still run on an updated dataset with the same columns.

Skill requirements:

- Include an ICA/open-question workflow:
  - start with a simple answerable question;
  - show data processing;
  - plot before modeling;
  - use methods already covered in course;
  - keep conclusions tied to the dataset;
  - avoid hard-coded absolute paths or identity-revealing paths.

## Review Findings: Corrections Needed in the User Template

The user-provided notes are valuable and should be used, but the skill should correct these points before turning them into reusable exam recipes:

- Use `pivot_longer()` and `pivot_wider()` as primary reshape functions; keep `gather()` only as a legacy equivalent.
- Do not blindly run `na.omit()`; first identify whether missingness is small and whether dropping rows changes the analysis.
- Do not remove duplicates automatically; duplicated patient/time records may be meaningful unless they are exact erroneous duplicates.
- For paired t-tests, check normality of paired differences.
- For Wilcoxon tests, use the correct name:
  - independent groups: Wilcoxon rank-sum / Mann-Whitney U;
  - paired groups: Wilcoxon signed-rank.
- For ANOVA power, Cohen's `f` must be `sqrt(eta_sq / (1 - eta_sq))`, not `eta_sq / (1 - eta_sq)`.
- For `pwr.anova.test`, `n` means per-group sample size; do not use `sum(x) / k` when `x` is the response values.
- For `pwr.t.test`, `n` for two-sample tests is per-group sample size; do not pass `n1 + n2` as if it were per group.
- For two-sample t-test effect size with unequal sample sizes, report that the simple pooled-SD Cohen's `d` is an approximation; Welch t-test may still be the better test.
- `chisq.test(cross_table$t)` is correct if `cross_table$t` is the count matrix; `chisq.test(cross_table)` is not.
- The chi-square assumptions refer to expected counts, not observed cell frequencies.
- Use "fail to reject H0" rather than "accept H0".
- Do not say "prove"; say "provide evidence", "consistent with", or "not enough evidence".
- The mock tooth-growth question is a two-factor ANOVA problem (`supp * dose`), so the skill must not limit ANOVA to one-way only.
- Optional packages such as `factoextra`, `CommKern`, `pwr`, `gmodels`, and `car` should have base-R fallbacks where possible.

## Constraints

- Do not modify the source course folders.
- Do not implement the skill bundle until this SPEC is accepted or revised.
- Keep extracted content as concise derived notes, not wholesale document copies.
- Preserve exam-use practicality: the output should optimize for speed, reproducibility, and answer correctness under time pressure.
- Prefer deterministic local extraction tools:
  - `pdftotext` for text PDFs.
  - Direct reading for `.Rmd`, `.R`, `.csv`, and `.txt`.
  - `officecli` or a local PPTX parser for slide text extraction if needed.
- Avoid network dependency.
- Use the user's local course materials as the authority when they conflict with generic statistics advice.

## Proposed Output

Create one self-contained skill bundle:

- `.codex/skills/ads2-online-exam-solver/SKILL.md`
- `.codex/skills/ads2-online-exam-solver/references/exam-answer-workflow.md`
- `.codex/skills/ads2-online-exam-solver/references/data-cleaning-and-reshaping.md`
- `.codex/skills/ads2-online-exam-solver/references/ads2-method-map.md`
- `.codex/skills/ads2-online-exam-solver/references/r-exam-recipes.md`
- `.codex/skills/ads2-online-exam-solver/references/probability-bayes-markov.md`
- `.codex/skills/ads2-online-exam-solver/references/power-and-sample-size.md`
- `.codex/skills/ads2-online-exam-solver/references/simulation-and-bootstrap.md`
- `.codex/skills/ads2-online-exam-solver/references/clustering-and-feature-extraction.md`
- `.codex/skills/ads2-online-exam-solver/references/interpretation-templates.md`
- `.codex/skills/ads2-online-exam-solver/references/common-traps.md`
- `.codex/skills/ads2-online-exam-solver/references/source-index.md`

Optional, only if useful after extraction:

- `.codex/skills/ads2-online-exam-solver/scripts/quick_data_audit.R`
- `.codex/skills/ads2-online-exam-solver/scripts/model_checklist.R`

The skill should trigger on phrases such as:

- ADS exam
- ADS2 exam
- Applied Data Science problem
- online problem solve exam
- solve this statistics question in R
- interpret this regression/correlation/ANOVA/categorical data output
- choose the right statistical test for this ADS question

## Skill Design

### `SKILL.md`

Purpose:

- Tell the agent when to use the skill.
- Establish the expected exam-solving workflow.
- Point to the relevant reference files only when needed.

Core workflow:

1. Start from the exam report workflow, not the statistical test.
2. Import, inspect, clean, and reshape data only as needed.
3. Classify the question type.
4. Identify variables, response/explanatory roles, measurement scale, sample structure, and independence assumptions.
5. Select the method using `ads2-method-map.md`.
6. Generate or inspect R code using `r-exam-recipes.md`.
7. Check assumptions and edge cases.
8. Interpret results in exam language using `interpretation-templates.md`.
9. Cross-check against `common-traps.md`.
10. Produce a concise final answer with method, code/output summary, conclusion, assumptions, limitations, and next-step recommendation.

### `ads2-method-map.md`

Expected content:

- Decision table from problem statement to method.
- Variable-type logic:
  - numeric vs numeric
  - numeric response vs categorical predictor
  - categorical response vs categorical predictor
  - binary/multiclass outcomes
  - repeated/grouped observations if present
- Common ADS2 methods likely represented in the materials:
  - data cleaning and reshaping
  - Bayes theorem and conditional probability
  - Markov-chain simulation
  - power and sample-size analysis
  - correlation
  - simple linear regression
  - multiple linear regression
  - t-tests
  - Wilcoxon rank-sum and signed-rank tests
  - Kruskal-Wallis test
  - one-way and two-way/factorial ANOVA
  - chi-square tests
  - Fisher's exact test
  - categorical-data summaries
  - confidence intervals and hypothesis tests
  - k-means clustering
  - basic feature extraction if supported by course files
- For each method:
  - Use when
  - Do not use when
  - R function pattern
  - Key assumptions
  - Output fields to read
  - Final-answer phrasing

### `r-exam-recipes.md`

Expected content:

- Copy-ready but adaptable R code patterns.
- Data loading and inspection:
  - `read.csv`
  - `str`
  - `summary`
  - `table`
  - `is.na`
  - factor conversion
- Exploratory plots:
  - scatterplot
  - boxplot
  - histogram
  - barplot/mosaic-style categorical summaries if used in course materials
- Modeling:
  - `t.test`
  - `wilcox.test`
  - `kruskal.test`
  - `aov`
  - `TukeyHSD`
  - `cor`
  - `cor.test`
  - `lm`
  - `summary(lm_fit)`
  - `confint`
  - `predict`
  - `anova`
  - `chisq.test`
  - `fisher.test`
  - `kmeans`
  - `power.t.test`
  - `pwr` package functions when installed, with base/simulation fallback
- Diagnostics:
  - residual plots
  - normality checks
  - leverage/outliers if course materials support it
- Output interpretation examples.

### `interpretation-templates.md`

Expected content:

- Short exam-answer templates for:
  - data cleaning and preprocessing
  - hypothesis statement
  - p-value decision
  - Bayes/posterior probability interpretation
  - power analysis and sample-size recommendation
  - confidence interval interpretation
  - slope interpretation
  - intercept interpretation
  - correlation interpretation
  - coefficient significance
  - model fit and R-squared
  - categorical association
  - predicted value interpretation
  - limitations and assumptions
  - next experimental step
- Templates should use placeholders, not hard-coded answers.

### `common-traps.md`

Expected content:

- Misreading response vs explanatory variable.
- Treating categorical variables as numeric without justification.
- Reporting p-values without conclusion in context.
- Saying "prove" instead of "evidence suggests".
- Interpreting intercepts outside meaningful range.
- Confusing correlation with causation.
- Forgetting units.
- Ignoring missing values.
- Using a paired test when samples are independent, or vice versa.
- Checking normality on raw paired samples instead of paired differences.
- Treating ANOVA as only one-way when the design has two factors.
- Ignoring interactions in two-factor designs.
- Using `pwr` functions with total sample size when per-group sample size is required.
- Reporting chi-square expected-count rules as observed-count rules.
- Running `chisq.test()` on a display object instead of a matrix/table.
- Overclaiming from model fit or small samples.
- Selecting a method from variable names rather than data type and question wording.

### `source-index.md`

Expected content:

- Table of source files used.
- For each source:
  - file path
  - file type
  - topic tags
  - extraction status
  - relevance level
  - notes on what was distilled

This is important so the skill remains auditable and can be updated later.

## Implementation Plan

### Phase 1: Source Inventory

Actions:

1. Build a complete file inventory for both source folders. Completed in second-pass review.
2. Filter to ADS-relevant files. Completed in second-pass review.
3. Categorize files into exam papers, practicals, slides, datasets, and unrelated material. Completed in second-pass review.
4. Create an initial `source-index.md` during implementation from the reviewed source list.

Acceptance criteria:

- Every ADS-relevant source file is listed.
- Non-ADS material is explicitly excluded or marked low relevance.
- The index can explain why each source was used or skipped.

### Phase 2: Text Extraction

Actions:

1. Extract PDF text with `pdftotext` into a temporary workspace under `/private/tmp`. Completed for review under `/private/tmp/ads2_skill_review`.
2. Extract R and RMarkdown source directly. Completed for review under `/private/tmp/ads2_skill_review/r_text`.
3. Extract PPTX text only for slides that fill gaps not covered by PDFs/Rmd. Completed for readable slides under `/private/tmp/ads2_skill_review/pptx_text`.
4. Summarize dataset columns using local commands or R, without copying full datasets. Partially completed for key exam datasets; finish during source-index creation.

Acceptance criteria:

- Each relevant source has either extracted text, parsed code, or a documented reason for skipping.
- Dataset summaries include column names, dimensions, and apparent variable types.
- No source file is modified.

### Phase 3: Topic Distillation

Actions:

1. Identify repeated exam question patterns from mock/past papers and Rmd solutions. Completed in this SPEC's review findings.
2. Map each pattern to:
   - statistical method
   - R code pattern
   - assumptions
   - interpretation language
   - likely traps
3. Compare practical/lecture content with exam files to fill missing theory and code details. Completed at the planning level; implementation should convert this into reference files.
4. Convert the user's provided notes into polished, corrected exam recipes.

Acceptance criteria:

- Method map covers all recurring ADS/ADS2 problem types observed in the materials.
- Each method has at least one code pattern and one interpretation template.
- Course-specific wording or conventions are preserved when visible in the sources.

### Phase 4: Skill Bundle Draft

Actions:

1. Create `.codex/skills/ads2-online-exam-solver/SKILL.md`.
2. Create reference files listed in the proposed output.
3. Keep `SKILL.md` compact and move detailed reference material into `references/`.
4. Add optional scripts only if repeated data-audit or model-checking code is clearly useful.
5. Include bilingual Chinese/English explanation where it preserves the user's exam-prep notes, while keeping R code and final exam report phrasing in clear English unless asked otherwise.

Acceptance criteria:

- Skill body stays readable and points to references by task.
- References are practical during an online exam.
- No reference file is bloated with raw course-document dumps.

### Phase 5: Validation

Actions:

1. Test the skill manually against at least six representative prompts:
   - Bayes/turtles posterior probability prompt.
   - k-means neuron classification prompt.
   - tooth-growth two-factor ANOVA prompt.
   - genotype/Mendelian chi-square prompt.
   - paired before/after t-test or Wilcoxon prompt.
   - ICA/open trend-analysis prompt.
2. Use past/mock exam questions where feasible.
3. Check that final answers are concise, contextual, and do not overclaim.
4. Check that the skill corrects the known risky template errors listed above.

Acceptance criteria:

- The skill improves answer structure compared with a generic statistics answer.
- R code is syntactically plausible and aligned with course examples.
- Interpretations include context, units where available, and uncertainty language.

### Phase 6: Optional Installation/Portability

Actions:

1. Keep the repo-local skill under `.codex/skills/`.
2. If the user wants global reuse, copy or mirror it into a user-level Codex skill location after approval.
3. If the user wants a shareable archive, create a portable folder or zip after approval.

Acceptance criteria:

- The local skill works in this repo.
- Any global install or external copy is done only after explicit user approval.

## Verification Commands

Planned local checks:

```bash
find .codex/skills/ads2-online-exam-solver -maxdepth 3 -type f | sort
```

```bash
sed -n '1,220p' .codex/skills/ads2-online-exam-solver/SKILL.md
```

```bash
rg -n "TODO|TBD|lorem|placeholder" .codex/skills/ads2-online-exam-solver
```

If R scripts are added:

```bash
Rscript --vanilla .codex/skills/ads2-online-exam-solver/scripts/quick_data_audit.R --help
```

## Risks and Mitigations

- Risk: PDFs are scanned or poorly extracted.
  - Mitigation: mark extraction gaps in `source-index.md`; only use reliable text unless OCR is explicitly approved.
- Risk: The bundle becomes too large to load during an exam.
  - Mitigation: keep `SKILL.md` short and split references by task.
- Risk: Generic statistics conventions conflict with course expectations.
  - Mitigation: prioritize local course examples and note any uncertainty.
- Risk: Past papers include unrelated courses in the same folder.
  - Mitigation: filter by filename and content before distillation.
- Risk: The agent gives overconfident answers under exam time pressure.
  - Mitigation: include an explicit final-answer checklist and common traps file.

## Open Questions

1. Should the final skill be installed only repo-locally under `.codex/skills/`, or also copied to a global Codex skill directory for use outside this website repo?
2. Should the skill answer in English only, Chinese only, or bilingual exam-support style?
3. Should the skill include compact R helper scripts, or should it remain pure Markdown references?
4. Do you want the skill to include an "exam mode" final-answer format, for example: Method, R code, Key output, Conclusion, Assumptions?

## Next Step After Approval

After this SPEC is accepted, start implementation from the reviewed extraction notes. The first concrete implementation action will be to create `.codex/skills/ads2-online-exam-solver/`, then write `SKILL.md`, `source-index.md`, and the reference files from the second-pass findings above.
