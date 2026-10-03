"""Runtime transport controls; source-only until compilation is admitted."""
import json
import subprocess
import sys
from decimal import Decimal


def reject_constant(value):
    raise AssertionError("bare nonfinite JSON token: " + value)


def require(value, message):
    if not value:
        raise AssertionError(message)


def normalization(value):
    require(isinstance(value, dict), "normalization object")
    flags = value["availability"]
    require(len(flags) == 6 and all(isinstance(x, bool) for x in flags), "explicit normalization availability")
    for key in ("q_min", "postcast_profile_q", "logdetH", "logdetC", "stationarity", "variance_sensitivity", "mean", "variance"):
        require((value[key] is not None) == flags[0], key + " completion scope")
    require((value["log_relative_integral"] is not None) == flags[3], "integral availability")
    require(value["cdf_nodes"] == value["cdf_attempts"] == value["bisections"] == 0, "no CDF work")


def main():
    binary = sys.argv[1]
    identity = ["0" * 40, "1" * 64, "2" * 64, "3" * 64]
    for mode in ("accepted", "schur-refused", "shape-refused", "wire-nonfinite-control", "wire-partial-layout-control"):
        run = subprocess.run([binary, *identity, mode], capture_output=True, text=True, check=False)
        require(run.returncode == 0, mode + " caller refusal/acceptance control")
        result = json.loads(run.stdout, parse_float=Decimal, parse_constant=reject_constant)
        if mode == "wire-nonfinite-control":
            require(result["scope"] == "serializer-only-synthetic-witness", "nonfinite control scope")
            require(result["witness"] == {"kind": "nonfinite", "classification": "positive-infinity"}, "typed actual synthetic witness")
            require(result["absent"] is None, "absent witness is null")
            continue
        require(result["schema_version"] == 1 and result["mode"] == mode, "schema identity")
        declared = result["declared_wrapper_identity"]
        require(list(declared.values()) == identity, "single-quoted owning identity fields")
        require(isinstance(result["preparation"]["availability"], list), "preparation availability")
        require(isinstance(result["output_layout_available"], bool), "outer output layout availability")
        if result["output_layout_available"]:
            count = len(result["training_normalizations"])
            require(all(len(result[key]) == count for key in (
                "training_refusal_witnesses", "training_offset_subtraction_rounding_estimates",
                "training_cross_projection_error_estimates", "training_offset_diagnostics_available",
                "training_cross_projection_diagnostics_available")), "earned matching training layout")
        if mode == "wire-partial-layout-control":
            require(result["batch_origin"] == "serializer-only-partial-allocation-layout", "partial layout control scope")
            require(not result["output_layout_available"], "partial layout unavailable")
            require(result["training_offset_subtraction_rounding_estimates"] == [None] and
                    result["training_cross_projection_error_estimates"] == [None], "missing flags withhold actual values")
            require(result["training_offset_diagnostics_available"] == [] and
                    result["training_cross_projection_diagnostics_available"] == [], "missing flags retained")
            continue
        for value in result["training_normalizations"]:
            normalization(value)
        for value in result["densities"]:
            require(isinstance(value["density_available"], bool), "density flag")
            require((value["log_density"] is not None) == value["density_available"], "reported density availability")
            normalization(value["joint_normalization"])
        if mode == "accepted":
            require(result["actual_source_order"] == ["r0", "r1", "r2"], "source order")
            require(result["actual_training_order"] == ["r0", "r1"], "training order")
            metadata = result["actual_metadata"]
            require(metadata["active_original_axes"] == [0, 2], "active original axis map")
            require(metadata["original_parameter_order"] == ["original0", "original1", "original2"], "point mass excluded from active order")
            require(metadata["fixed_original_axis"] == 1 and metadata["fixed_value"] == 0, "point mass witness")
            require(result["densities"][0]["density_available"], "accepted density")
        elif mode == "schur-refused":
            require(result["preparation"]["status"] == 4 and result["source_after_prepare_status"] == 0, "numerical refusal preserves source")
            require(result["densities"] == [] and result["preparation"]["refusal_witness"] is not None, "failed setup earned witness retained")
        else:
            require(result["batch_status"] == 2 and result["densities"] == [], "invalid pool withheld")


if __name__ == "__main__":
    main()
