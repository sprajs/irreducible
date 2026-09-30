"""Current structural request schema; scientific admission remains native."""
from copy import deepcopy


def build_schema(abi, fixture):
    fixture = deepcopy(fixture)
    fixture["properties"]["schema_version"] = {"const": 2}
    fixture["$id"] = "irreducible.fixture-run.v2"
    number = {"type": "number"}
    text = {"type": "string"}
    uint = {"type": "integer", "minimum": 0, "maximum": 18446744073709551615}

    def obj(properties, optional=()):
        return {"type": "object", "additionalProperties": False,
                "properties": properties,
                "required": [key for key in properties if key not in optional]}

    def array(items, maximum=None):
        result = {"type": "array", "items": items}
        if maximum is not None:
            result["maxItems"] = maximum
        return result

    def enum(values):
        return {"enum": list(values)}

    def outputs(values):
        return {**array(enum(values)), "minItems": 1, "uniqueItems": True}

    def ref(name):
        return {"$ref": "#/$defs/" + name}

    def nullable(value):
        return {"anyOf": [value, {"type": "null"}]}

    def operation(name, fields, optional=(), version=2):
        return obj({"schema_version": {"const": version},
                    "operation": {"const": name}, **fields}, optional)

    def variant(kind, fields=None):
        return obj({"kind": {"const": kind}, **(fields or {})})

    arithmetic = enum(["binary64", "wide"])
    geometry = variant("flat_flrw")
    expansion = {"oneOf": [variant("lcdm", {"omega_m": number}),
        variant("constant_q", {"q": number}),
        variant("cpl", {"omega_m": number, "w0": number, "wa": number}),
        variant("fixed_q5", {"q": {**array(number), "minItems": 5, "maxItems": 5}})]}
    effect = {"oneOf": [variant("none"), variant("grey_log1p_magnitude", {"epsilon_mag": number})]}
    observer = obj({"redshift": number, "convention": enum([
        "geometric_same_redshift", "released_zhd_zhel"])})
    integration = obj({"absolute_tolerance": {**number, "minimum": 0},
        "relative_tolerance": {**number, "minimum": 0},
        "maximum_evaluations": uint,
        "maximum_depth": {"type": "integer", "minimum": 0, "maximum": 60}})
    projection = obj({"integration": nullable(integration), "maximum_queries": uint,
        "maximum_callbacks": uint, "maximum_segment_visits": uint}, ("integration",))
    preparation = {"arithmetic": arithmetic, "maximum_matrix_elements": uint,
        "maximum_string_bytes": uint, "maximum_native_bytes": uint,
        "maximum_forward_sensitivity": {**number, "exclusiveMinimum": 0}}
    evaluation = obj({"arithmetic": arithmetic, "projection": projection,
        "maximum_models": uint, "maximum_array_elements": uint,
        "maximum_native_bytes": uint,
        "maximum_forward_sensitivity": {**number, "exclusiveMinimum": 0}})
    tags = abi["observation_tags"]
    metadata = obj({key: enum(values) for key, values in tags.items()
                    if key not in ("selection", "status")})
    observations = operation("observations.prepare", {
        "table": text, "uncertainty": nullable(text), "metadata": ref("observation_metadata"),
        "resources": obj({"maximum_asset_bytes": {**uint, "maximum": 268435456},
            "maximum_rows": {**uint, "maximum": 10000},
            "maximum_matrix_elements": {**uint, "maximum": 25000000},
            "maximum_string_bytes": {**uint, "maximum": 16777216},
            "maximum_preparation_bytes": uint}),
        "exports": {**array(enum(["source_values", "covariance", "original_fields"])), "uniqueItems": True},
        **{key: text for key in ("calibration_provenance", "dependence_provenance",
                                "quality_dictionary", "ordering_provenance")},
        "source_selection": array(enum([0, 1]), 10000)}, ("uncertainty",))
    source_common = {key: text for key in ("ordering_provenance", "calibration_provenance",
        "dependence_provenance", "redshift_convention", "ruler_convention")}
    bao_source = {"oneOf": [obj({"profile": {"const": "desi_dr2_all_gccomb_13_v1"},
        "mean": text, "covariance": text, "maximum_asset_bytes": uint, **source_common}),
        obj({"profile": {"const": "synthetic_inline_v1"},
            "rows": array(obj({"id": text, "z": number,
                "observable": enum(["DM_over_rs", "DH_over_rs", "DV_over_rs"]), "value": number}), 4096),
            "covariance": array(number, 16777216),
            "covariance_axis_ids": array(text, 4096), **source_common})]}
    selection = {"oneOf": [variant("pantheon_zhd_gt001"), variant("explicit", {
        "source_indices": array(uint, 4096), "coordinates": array(obj({
            "z_expansion": number, "observer": observer}), 4096)})]}
    defs = {"fixture": fixture, "observation_metadata": metadata,
        "observations": observations, "expansion_spec": expansion,
        "geometry": geometry, "source_effect": effect, "bao_source": bao_source}
    defs["background"] = operation("background.evaluate", {
        "geometry": geometry, "models": array(expansion, 65536),
        "queries": array(obj({"z_expansion": number,
            "requested_outputs": outputs(["expansion", "radial", "luminosity_shape", "clock", "physical", "kinematics"]),
            "observer": nullable(observer), "physical_scale": nullable(obj({"h0_km_s_mpc": number}))},
            ("observer", "physical_scale")), 65536),
        "numerical_policy": obj({"integration": nullable(integration),
            **{key: uint for key in ("maximum_models", "maximum_queries", "maximum_slots",
                 "maximum_callbacks", "maximum_segment_visits", "maximum_native_bytes")}}, ("integration",))})
    defs["supernova"] = operation("supernova.profile", {
        "observations": ref("observations"), "selection": selection,
        "preparation_policy": obj({**preparation, "maximum_selected_rows": uint}),
        "models": array(obj({"expansion": expansion, "geometry": geometry, "source_effect": effect}), 65536),
        "requested_outputs": outputs(["score", "geometric_shape", "magnitude_effect",
            "corrected_residuals", "profiled_residuals", "diagnostics"]),
        "numerical_policy": evaluation})
    defs["bao"] = operation("bao.density", {"source": bao_source,
        "preparation_policy": obj({**preparation, "maximum_queries": uint}),
        "models": array(obj({"expansion": expansion, "geometry": geometry, "h0_rd_km_s": number}), 65536),
        "requested_outputs": outputs(["normalized_density", "predictions", "residuals"]),
        "numerical_policy": evaluation})
    prior = obj({"response": array(number), "ordered_ids": array(text), "mean": number,
        "variance": number, "latent_identity": text, "independence_declared": {"type": "boolean"}})
    defs["statistics"] = operation("statistics.gaussian", {
        "observations": ref("observations"), "mode": enum(["normalized_density", "profile_offset_score"]),
        "residual_unit": text, "response_unit": text, "ordered_ids": array(text),
        "residuals": array(array(number)), "response": nullable(array(number)),
        "proper_prior": nullable(prior), "selection": enum(tags["selection"]),
        "maximum_preparation_bytes": uint, "maximum_evaluation_bytes": uint,
        "maximum_batch_elements": uint,
        "maximum_forward_sensitivity": {**number, "exclusiveMinimum": 0}}, ("response", "proper_prior"))
    defs["quantity_metadata"] = obj({key: enum(values) for key, values in abi["quantity_tags"].items()
                                     if key != "quantity_status"})
    defs["quantity"] = operation("quantity.convert", {"values": array(number),
        "source": ref("quantity_metadata"), "target": ref("quantity_metadata")})
    defs["numerics"] = operation("numerics.scalar_batch", {
        "method": enum(abi["numerical_tags"]["operation"]),
        "values": array(number, abi["max_batch_elements"])})
    defs["sound_horizon"] = operation("cosmology.sound_horizon", {
        "points": array(obj({"h0_km_s_mpc": number, "omega_m": number,
            "omega_r": number, "omega_b": number, "omega_gamma": number,
            "z_drag": number, "drag_origin": text}), 65536),
        "numerical_policy": obj({"absolute_tolerance_mpc": number,
            "relative_tolerance": number, "maximum_callbacks_per_point": uint,
            "maximum_depth": {**uint, "maximum": 60}, "maximum_points": uint,
            "maximum_total_callbacks": uint, "maximum_native_bytes": {**uint, "maximum": 1073741824}})})
    return {"$schema": "https://json-schema.org/draft/2020-12/schema",
        "$id": "cosmology.run.v2", "$defs": defs,
        "oneOf": [ref(name) for name in ("fixture", "quantity", "numerics", "observations",
                                          "background", "statistics", "supernova", "bao", "sound_horizon")],
        "description": "Structural compiled requests; native domains, numerical gates and scientific qualifications remain distinct."}
