# Select fold-ranked CpGs and fit survival models.

import argparse
import multiprocessing as mp
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import warnings

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler


METHODS = ("beta", "chalm", "camda")
TOP_NS = (5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)
RIDGE_PENALTY = 0.1
PROJECTS = (
    "TCGA-BLCA", "TCGA-BRCA", "TCGA-COAD", "TCGA-GBM", "TCGA-KICH",
    "TCGA-KIRC", "TCGA-KIRP", "TCGA-LGG", "TCGA-LIHC", "TCGA-LUAD",
    "TCGA-LUSC", "TCGA-PAAD", "TCGA-PRAD", "TCGA-READ", "TCGA-SKCM",
)
MODELS = (("beta",), ("chalm", "camda"), METHODS)
CONTEXT = None


def phenotype(path, project):
    data = pd.read_csv(path, low_memory=False)
    data = data.loc[data["GSE_ID"].eq(project)].copy()
    data["time"] = pd.to_numeric(data["Overall_Survival(Year)"], errors="coerce")
    vital = data["Vital_Status"].astype("string").str.strip().str.lower()
    data["status"] = vital.map({"dead": 1, "alive": 0})
    data = data.dropna(subset=["time", "status"])
    data = data.loc[data["time"].gt(0)]
    result = data.set_index("GSM_ID")[["time", "status"]]
    if not result.index.is_unique:
        raise ValueError("Duplicated sample IDs in phenotype")
    return result


def cox_ranking(path):
    table = pd.read_csv(path, usecols=["cpg", "pValue"], low_memory=False)
    if table["cpg"].isna().any() or table["cpg"].duplicated().any():
        raise ValueError("Missing or duplicated CpG IDs in fold Cox output")
    table["pValue"] = pd.to_numeric(table["pValue"], errors="coerce")
    table = table.loc[np.isfinite(table["pValue"]) & table["pValue"].between(0, 1)]
    table["cpg"] = table["cpg"].astype(str)
    table = table.sort_values("pValue", kind="mergesort").reset_index(drop=True)
    return table["cpg"].tolist()


def matrix_path(base_dir, camda_dir, method, project):
    if method == "beta":
        return base_dir / "beta" / (project + "_assay.tsv")
    if method == "chalm":
        path = base_dir / "chalm_camda" / (project + ".chalm.filtered.tsv")
        fallback = base_dir / "chalm_camda" / (project + ".chalm.tsv")
        return path if path.is_file() else fallback
    paths = (
        camda_dir / (project + ".camda.filtered.tsv"),
        camda_dir / (project + ".camda.tsv"),
        camda_dir / project / (project + ".camda.tsv"),
    )
    return next((path for path in paths if path.is_file()), paths[0])


def load_matrix(path, selected):
    if not path.is_file():
        return None
    table = pd.read_csv(path, sep="\t", index_col=0, low_memory=False)
    table.index = table.index.astype(str)
    table.columns = [str(value).replace(".", "-") for value in table.columns]
    if not table.columns.is_unique:
        raise ValueError("Duplicated sample IDs in prediction matrix")
    common = table.index.intersection(list(selected))
    return table.loc[common].T


def evaluate_fold(x, pheno, train_idx, test_idx):
    x_train, x_test = x.iloc[train_idx], x.iloc[test_idx]
    t_train, t_test = pheno.time.iloc[train_idx], pheno.time.iloc[test_idx]
    e_train, e_test = pheno.status.iloc[train_idx], pheno.status.iloc[test_idx]
    row = {
        "status": "fit_failed", "c_index": np.nan, "error": "",
        "n_features": x.shape[1], "penalty": np.nan,
        "n_train": len(train_idx), "n_test": len(test_idx),
        "test_events": int(e_test.sum()),
    }
    try:
        imputer = SimpleImputer(strategy="mean")
        scaler = StandardScaler()
        x_train = scaler.fit_transform(imputer.fit_transform(x_train))
        x_test = scaler.transform(imputer.transform(x_test))
        train = pd.DataFrame(x_train, columns=["f{}".format(i) for i in range(x_train.shape[1])])
        train["T"], train["E"] = t_train.to_numpy(), e_train.to_numpy()
        penalty = RIDGE_PENALTY
        model = CoxPHFitter(penalizer=penalty)
        model.fit(train, duration_col="T", event_col="E")
        test = pd.DataFrame(x_test, columns=train.columns[:-2])
        risk = model.predict_partial_hazard(test)
        row.update(
            status="ok",
            c_index=float(concordance_index(t_test, -risk, e_test)),
            n_features=int(x_train.shape[1]),
            penalty=float(penalty),
        )
    except Exception as error:
        row["error"] = "{}: {}".format(type(error).__name__, error)
    return row


def predict_model(task):
    method_set, top_n = task
    state = CONTEXT
    project = state["project"]
    pheno = state["pheno"]
    roster = state["roster"]
    assignment = state["assignment"]
    rankings = state["rankings"]
    frames = state["frames"]
    method_name = "+".join(method_set)
    rows = []
    repeat_scores = []
    sample_index = {sample: index for index, sample in enumerate(roster)}

    for repeat in range(1, 4):
        fold_scores = []
        for fold in range(1, 6):
            scope = "r{}_f{}".format(repeat, fold)
            actual = assignment.loc[assignment["scope"].eq(scope)]
            if actual["Sample"].tolist() != roster:
                raise ValueError("Fold sample order changed")
            train_idx = np.array([
                sample_index[sample] for sample in actual.loc[
                    actual["partition"].eq("train"), "Sample"
                ]
            ])
            test_idx = np.array([
                sample_index[sample] for sample in actual.loc[
                    actual["partition"].eq("test"), "Sample"
                ]
            ])

            parts = []
            for method in method_set:
                frame = frames[method]
                chosen = [
                    cpg for cpg in rankings[scope][method][:top_n]
                    if frame is not None and cpg in frame.columns
                ]
                if len(chosen) != top_n:
                    raise ValueError(
                        "{} {} Top-{} features are incomplete in {}".format(
                            project, method, top_n, scope
                        )
                    )
                part = frame[chosen].copy()
                if len(method_set) > 1:
                    part.columns = ["{}_{}".format(cpg, method) for cpg in part.columns]
                parts.append(part)

            result = {
                "Project": project, "Method": method_name, "TopN": top_n,
                "repeat_id": repeat, "outer_fold": fold, "scope": scope,
                "status": "skipped_no_features", "c_index": np.nan, "error": "",
                "n_features": 0, "penalty": np.nan,
                "n_train": len(train_idx), "n_test": len(test_idx),
                "test_events": int(pheno.status.iloc[test_idx].sum()),
            }
            if not parts:
                rows.append(result)
                continue
            x = pd.concat(parts, axis=1) if len(parts) > 1 else parts[0]
            x = x.reindex(roster)
            if len(parts) > 1:
                x = x.dropna(axis=1, how="all")
                if len(x) <= 20:
                    result["status"] = "skipped_small_combination"
                    rows.append(result)
                    continue
            if len(x) < 10 or x.shape[1] == 0:
                result["status"] = "skipped_precheck"
                rows.append(result)
                continue
            fit = evaluate_fold(x, pheno, train_idx, test_idx)
            result.update(fit)
            rows.append(result)
            if result["status"] == "ok":
                fold_scores.append(result["c_index"])
        if fold_scores:
            repeat_scores.append(float(np.mean(fold_scores)))

    failed = [row for row in rows if row["status"] != "ok"]
    if len(rows) != 15 or failed:
        raise RuntimeError(
            "{} {} Top-{} did not complete all 15 folds: {}".format(
                project, method_name, top_n,
                "; ".join(
                    "{}={}".format(row["scope"], row["status"])
                    for row in failed
                )
            )
        )
    c_index = float(np.mean([row["c_index"] for row in rows]))
    summary = {
        "Project": project, "Method": method_name, "TopN": top_n,
        "C_index": c_index, "valid_repeats": len(repeat_scores),
        "valid_folds": sum(row["status"] == "ok" for row in rows),
    }
    return summary, rows


def fit_project(project, phenotype_path, base_data_dir, camda_data_dir,
                fold_assignments, cox_dir, out_dir, threads):
    args = argparse.Namespace(
        project=project, phenotype=phenotype_path,
        base_data_dir=base_data_dir, camda_data_dir=camda_data_dir,
        fold_assignments=fold_assignments, cox_dir=cox_dir,
        out_dir=out_dir, threads=threads,
    )
    if not 1 <= args.threads <= 5:
        raise ValueError("Use between one and five threads")

    assignment = pd.read_csv(args.fold_assignments, dtype={"Sample": str})
    scopes = {"r{}_f{}".format(r, f) for r in range(1, 4) for f in range(1, 6)}
    if set(assignment["scope"]) != scopes or assignment.duplicated(["scope", "Sample"]).any():
        raise ValueError("Expected 3 repeats x 5 folds")
    roster = assignment.loc[assignment["scope"].eq("r1_f1"), "Sample"].tolist()
    for scope in scopes:
        if assignment.loc[assignment["scope"].eq(scope), "Sample"].tolist() != roster:
            raise ValueError("Sample order differs across folds")

    rankings = {}
    for scope in sorted(scopes):
        rankings[scope] = {}
        for method in METHODS:
            path = args.cox_dir / args.project / scope / (
                "{}_{}_cox_purity_immune.csv".format(args.project, method)
            )
            rankings[scope][method] = cox_ranking(path)

    pheno = phenotype(args.phenotype, args.project)
    if not set(roster).issubset(pheno.index):
        raise ValueError("Phenotype does not contain every fold sample")
    pheno = pheno.loc[roster]

    frames = {}
    for method in METHODS:
        selected = set()
        for scope in scopes:
            selected.update(rankings[scope][method][:max(TOP_NS)])
        frame = load_matrix(matrix_path(args.base_data_dir, args.camda_data_dir,
                                        method, args.project), selected)
        if frame is not None:
            frame = frame.loc[frame.index.intersection(pheno.index)]
            if frame.index.tolist() != roster:
                raise ValueError("Prediction matrix sample order differs for " + method)
        frames[method] = frame

    global CONTEXT
    CONTEXT = dict(project=args.project, pheno=pheno, roster=roster,
                   assignment=assignment, rankings=rankings, frames=frames)
    tasks = [(model, n) for model in MODELS for n in TOP_NS]
    if args.threads == 1:
        results = list(map(predict_model, tasks))
    else:
        with mp.get_context("fork").Pool(args.threads) as pool:
            results = pool.map(predict_model, tasks)
    CONTEXT = None

    summary = pd.DataFrame([item[0] for item in results])
    folds = pd.DataFrame([row for item in results for row in item[1]])
    if len(summary) != 33 or len(folds) != 495:
        raise ValueError("Incomplete three-model, 11-TopN, 15-fold grid")
    args.out_dir.mkdir(parents=True, exist_ok=True)
    summary.to_csv(args.out_dir / "model_coverage.csv", index=False)
    folds.to_csv(args.out_dir / "fold_results.csv", index=False)
    return summary


def run_cox(project, threads, paths, project_root, log_path):
    rscript = os.environ.get("RSCRIPT") or shutil.which("Rscript")
    if not rscript:
        raise ValueError("Set RSCRIPT or make Rscript available on PATH")
    code = Path(__file__).resolve().parent
    command = [
        rscript, str(code / "tcga_fold_cox.R"),
        "--projects", project, "--cores", str(threads),
        "--phenotype_dir", str(paths["gdc_pheno"]),
        "--survival_phenotype", str(paths["phenotype"]),
        "--beta_dir", str(paths["base"] / "beta"),
        "--chalm_dir", str(paths["base"] / "chalm_camda"),
        "--camda_dir", str(paths["camda"]),
        "--ref_dir", str(paths["references"]),
        "--outdir", str(project_root),
    ]
    with log_path.open("w", encoding="utf-8") as log:
        subprocess.run(command, check=True, stdout=log, stderr=subprocess.STDOUT)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--project", choices=PROJECTS)
    group.add_argument("--all", action="store_true")
    parser.add_argument("--threads", type=int, default=5)
    args = parser.parse_args()
    if not 1 <= args.threads <= 5:
        raise ValueError("Use between one and five threads")

    paths = {
        "phenotype": Path(os.environ["SURVIVAL_PHENOTYPE"]),
        "gdc_pheno": Path(os.environ["GDC_PHENOTYPE_DIR"]),
        "base": Path(os.environ["BASE_DATA_DIR"]),
        "camda": Path(os.environ["CAMDA_DATA_DIR"]),
        "references": Path(os.environ["EPISCORE_REF_DIR"]),
    }
    output_root = Path(tempfile.gettempdir()) / "tcga_survival"
    projects = PROJECTS if args.all else (args.project,)
    tables = []
    for project in projects:
        project_root = output_root / "projects" / project
        folds = project_root / "folds" / "split_assignments.csv"
        cox = project_root / "cox"
        prediction = project_root / "prediction"
        logs = project_root / "logs"
        logs.mkdir(parents=True, exist_ok=True)
        run_cox(project, args.threads, paths, project_root, logs / "cox.log")
        table = fit_project(
            project, paths["phenotype"], paths["base"], paths["camda"],
            folds, cox, prediction, args.threads,
        )
        tables.append(table)
        print("[COMPLETE]", project)
    if args.all:
        all_results = pd.concat(tables, ignore_index=True)
        if set(all_results["Project"]) != set(PROJECTS):
            raise ValueError("Pan-cancer summary requires all 15 TCGA projects")
        all_results.to_csv(output_root / "model_coverage.csv", index=False)
        pan_cancer = all_results.groupby(["Method", "TopN"]).agg(
            median=("C_index", "median"),
            q25=("C_index", lambda values: values.quantile(0.25)),
            q75=("C_index", lambda values: values.quantile(0.75)),
            n_cancers=("Project", "nunique"),
        ).reset_index()
        if not pan_cancer["n_cancers"].eq(15).all():
            raise ValueError("Each pan-cancer estimate must include 15 cancer types")
        pan_cancer.to_csv(output_root / "pan_cancer_summary.csv", index=False)


if __name__ == "__main__":
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        main()
