import json
import os
from copy import deepcopy
from datetime import datetime


RUNS_DIR = os.path.join("results", "runs")
GENERATIONS_FILENAME = "generations.json"


from analysis.diversity import compute_population_diversity


def summarize_crossover_stats(stats):
    attempts = stats.get("attempts", 0)
    successes = stats.get("successes", 0)
    failures = stats.get("failures", 0)
    if attempts == 0:
        return {
            "success_rate_pct": None,
            "failure_rate_pct": None,
        }
    return {
        "success_rate_pct": round(100.0 * successes / attempts, 1),
        "failure_rate_pct": round(100.0 * failures / attempts, 1),
    }


def summarize_mutation_stats(stats):
    invocations = stats.get("invocations", 0)
    successes = stats.get("successes", 0)
    if invocations == 0:
        return {
            "success_rate_pct": None,
            "fallback_rate_pct": None,
            "avg_attempts_on_success": None,
        }
    return {
        "success_rate_pct": round(100.0 * successes / invocations, 1),
        "fallback_rate_pct": round(100.0 * stats.get("fallbacks", 0) / invocations, 1),
        "avg_attempts_on_success": round(
            stats.get("inner_attempts_on_success", 0) / successes,
            1,
        ) if successes else None,
    }


class EvolutionStatistics:
    """Stores metrics collected during a genetic algorithm run and exports JSON."""

    def __init__(self):
        self.reset()

    def reset(self):
        self.run_id = None
        self.started_at = None
        self.config = {}
        self.generations_data = []
        self.operators_summary = None
        self._next_generation = 1

    @property
    def generation_indices(self):
        return [entry["generation"] for entry in self.generations_data]

    @property
    def best_fitness(self):
        return [entry["best_fitness"] for entry in self.generations_data]

    @property
    def average_fitness(self):
        return [entry["average_fitness"] for entry in self.generations_data]

    def has_data(self):
        return len(self.generations_data) > 0

    def final_generation_entry(self):
        if not self.generations_data:
            return None
        return self.generations_data[-1]

    def begin_run(self, config):
        self.reset()
        self.run_id = datetime.now().strftime("%Y%m%d_%H%M%S")
        self.started_at = datetime.now().isoformat(timespec="seconds")
        self.config = deepcopy(config)
        self._write_export()

    def record_generation(
        self,
        population,
        crossover_mode=None,
        mutation_mode=None,
        crossover_stats=None,
        mutation_stats=None,
    ):
        if not population:
            return

        fitness_values = [individual.getQED() for individual in population]
        crossover_raw = deepcopy(crossover_stats or {})
        mutation_raw = deepcopy(mutation_stats or {})
        entry = {
            "generation": self._next_generation,
            "best_fitness": max(fitness_values),
            "average_fitness": sum(fitness_values) / len(fitness_values),
            "population_size": len(population),
            "diversity": compute_population_diversity(population),
            "crossover": {
                "mode": crossover_mode.name if crossover_mode is not None else None,
                "stats": crossover_raw,
                "summary": summarize_crossover_stats(crossover_raw),
            },
            "mutation": {
                "mode": mutation_mode.name if mutation_mode is not None else None,
                "stats": mutation_raw,
                "summary": summarize_mutation_stats(mutation_raw),
            },
        }
        self.generations_data.append(entry)
        self._next_generation += 1
        self._write_export()

    def update_operators_summary(
        self,
        crossover_mode,
        mutation_mode,
        crossover_stats,
        mutation_stats,
    ):
        crossover_raw = deepcopy(crossover_stats or {})
        mutation_raw = deepcopy(mutation_stats or {})
        self.operators_summary = {
            "updated_at": datetime.now().isoformat(timespec="seconds"),
            "crossover": {
                "mode": crossover_mode.name if crossover_mode is not None else None,
                "stats": crossover_raw,
                "summary": summarize_crossover_stats(crossover_raw),
            },
            "mutation": {
                "mode": mutation_mode.name if mutation_mode is not None else None,
                "stats": mutation_raw,
                "summary": summarize_mutation_stats(mutation_raw),
            },
        }
        self._write_export()

    def export_dict(self):
        payload = {
            "run_id": self.run_id,
            "started_at": self.started_at,
            "config": self.config,
            "generations": deepcopy(self.generations_data),
        }
        if self.operators_summary is not None:
            payload["operators_summary"] = deepcopy(self.operators_summary)
        return payload

    def run_directory(self):
        if not self.run_id:
            return None
        return os.path.join(RUNS_DIR, self.run_id)

    def export_path(self):
        run_dir = self.run_directory()
        if run_dir is None:
            return None
        return os.path.join(run_dir, GENERATIONS_FILENAME)

    def _write_export(self):
        if not self.run_id:
            return
        run_dir = self.run_directory()
        os.makedirs(run_dir, exist_ok=True)
        with open(self.export_path(), "w", encoding="utf-8") as export_file:
            json.dump(self.export_dict(), export_file, indent=2)

    @classmethod
    def from_export(cls, path):
        with open(path, encoding="utf-8") as export_file:
            data = json.load(export_file)

        stats = cls()
        stats.run_id = data.get("run_id")
        stats.started_at = data.get("started_at")
        stats.config = data.get("config", {})
        stats.generations_data = data.get("generations", [])
        stats.operators_summary = data.get("operators_summary")
        if stats.generations_data:
            stats._next_generation = stats.generations_data[-1]["generation"] + 1
        return stats
