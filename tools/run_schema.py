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
    thermal_text = {"type": "string", "minLength": 1, "maxLength": 256}
    thermal_species = obj({key: number for key in ("mass_ev", "temperature_today_kelvin", "statistical_weight")})
    thermal_physical = obj({**{key: number for key in ("h0_km_s_mpc", "physical_baryon_density",
        "physical_cdm_density", "tcmb_kelvin", "physical_massless_nonphoton_density")},
        "species": array(thermal_species, 16)})
    thermal_momentum = obj({"absolute_tolerance": number, "relative_tolerance": number,
        "maximum_callbacks_per_evaluation": {**uint, "maximum": 200000},
        "maximum_total_callbacks": {**uint, "maximum": 100000000},
        "maximum_depth": {**uint, "maximum": 64},
        "maximum_native_bytes": {**uint, "maximum": 67108864},
        "momentum_method": {"const": "direct_adaptive"}})
    thermal_prediction = obj({**{key: number for key in ("absolute_tolerance_mpc", "relative_tolerance",
        "absolute_tolerance_ratio", "relative_tolerance_ratio")},
        "maximum_callbacks_per_point": {**uint, "maximum": 20000000},
        "maximum_total_callbacks": {**uint, "maximum": 100000000},
        "maximum_depth": {**uint, "maximum": 64},
        "maximum_native_bytes": {**uint, "maximum": 67108864}, "thermal": thermal_momentum})
    thermal_caps = {"maximum_rows":64,"maximum_matrix_elements":4096,"maximum_models":16,
        "maximum_species_per_model":16,"maximum_string_bytes":65536,"maximum_output_array_elements":2048,
        "maximum_native_bytes":67108864,"maximum_preparation_native_bytes":67108864,
        "maximum_evaluation_native_bytes":67108864,"maximum_total_callbacks":200000000}
    defs["bao_thermal"] = operation("bao.thermal_density", {
        "source_semantics": {"const":"synthetic_controls"},
        "observation": obj({"ordered_row_ids":array(thermal_text,64),
            "queries":array(obj({"redshift":number,"observable":enum(["DM_over_rs","DH_over_rs","DV_over_rs"])}),64),
            "observed_ratios":array(number,64),"covariance_axis_ids":array(thermal_text,64),
            "covariance_row_major":array(number,4096),"ratio_unit":{"const":"one"},
            "covariance_unit":{"const":"ratio_squared"},**{key:thermal_text for key in (
                "table_identity","covariance_identity","ordering_provenance","calibration_provenance",
                "dependence_provenance","redshift_convention")}}),
        "models":array(obj({"id":thermal_text,"physical_model":thermal_physical,"z_drag":number,
            "drag_origin":thermal_text,"source_origin":thermal_text}),16),
        "requested_outputs":outputs(["normalized_density","predictions","residuals"]),
        "resource_policy":obj({"arithmetic":arithmetic,**{key:{**uint,"maximum":maximum} for key,maximum in thermal_caps.items()},
            "maximum_forward_sensitivity":number,"maximum_projection_log_density_error":number,"predictions":thermal_prediction})})
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
    posterior_text = {"type": "string", "minLength": 1, "maxLength": 256,
                      "description": "1..256 UTF-8 bytes; native byte admission also applies"}
    posterior_ids = array(posterior_text, 65536)
    posterior_numbers = array(number, 1000000)
    defs["gaussian_posterior"] = operation("statistics.gaussian_posterior", {
        "source_semantics": {"const": "synthetic_controls"},
        "noise": obj({"ordered_row_ids": posterior_ids, "event_ids": posterior_ids,
            "residual_unit": posterior_text, "covariance_row_major": posterior_numbers,
            "noise_identity": posterior_text, "calibration_identity": posterior_text,
            "dependence_identity": posterior_text, "ordering_provenance": posterior_text}),
        "design": obj({"ordered_row_ids": posterior_ids, "ordered_parameter_ids": posterior_ids,
            "parameter_units": posterior_ids, "column_units": posterior_ids,
            "values_row_major": posterior_numbers, "design_identity": posterior_text}),
        "parameter_prior": obj({"ordered_parameter_ids": posterior_ids, "parameter_units": posterior_ids,
            "shared_nuisance_ids": posterior_ids, "mean": posterior_numbers,
            "covariance_row_major": posterior_numbers, "prior_identity": posterior_text,
            "parameter_measure": posterior_text, "dependence_identity": posterior_text,
            "noise_independence_declared": {"type": "boolean"}}),
        "conditioning": obj({"ordered_row_ids": posterior_ids, "event_ids": posterior_ids,
            "case_ids": posterior_ids, "vectors": array(posterior_numbers, 65536)}),
        "resource_policy": obj({"arithmetic": arithmetic,
            "maximum_elements": {**uint, "maximum": 1000000},
            "maximum_cases": {**uint, "maximum": 65536},
            "maximum_string_bytes": {**uint, "maximum": 1048576},
            "maximum_native_bytes": {**uint, "maximum": 268435456},
            "maximum_work_units": {**uint, "maximum": 100000000},
            "maximum_forward_sensitivity": {**number, "exclusiveMinimum": 0}})})
    # Same bounded structural types; the predictive operation adds object/mask gates.
    posterior = defs["gaussian_posterior"]["properties"]
    predictive = operation("statistics.gaussian_predictive", {
        "source_semantics": {"const": "synthetic_controls"},
        "noise": posterior["noise"], "design": posterior["design"],
        "parameter_prior": posterior["parameter_prior"],
        "conditioning": posterior["conditioning"],
        "future_noise": posterior["noise"], "future_response": posterior["design"],
        "prediction": obj({"conditioning_identity": posterior_text,
            "dependence_identity": posterior_text, "future_covariance_unit": posterior_text,
            "future_measure": posterior_text,
            "future_noise_independence_declared": {"type": "boolean"},
            "noise_conditional_on_parameters_declared": {"type": "boolean"}}),
        "outputs": obj({"means": {"type": "boolean"}, "joint_log_densities": {"type": "boolean"}}),
        "future_vectors": obj({"ordered_row_ids": posterior_ids, "event_ids": posterior_ids,
            "vectors": array(posterior_numbers, 65536)}),
        "resource_policy": posterior["resource_policy"]}, ("future_vectors",))
    predictive["properties"]["resource_policy"] = {**posterior["resource_policy"],
        "properties": {**posterior["resource_policy"]["properties"],
            "maximum_forward_sensitivity": {**number, "exclusiveMinimum": 0, "maximum": 1e-10}}}
    predictive["allOf"] = [
        {"anyOf": [{"properties": {"outputs": {"properties": {"means": {"const": True}}}}},
                   {"properties": {"outputs": {"properties": {"joint_log_densities": {"const": True}}}}}]},
        {"if": {"properties": {"outputs": {"properties": {"joint_log_densities": {"const": True}}}}},
         "then": {"required": ["future_vectors"]}, "else": {"not": {"required": ["future_vectors"]}}}]
    defs["gaussian_predictive"] = predictive
    photometry_fields = ("luminosity_watt_per_metre", "rest_lower_metre", "rest_upper_metre",
        "observed_lower_metre", "observed_upper_metre", "luminosity_distance_metre", "redshift",
        "collecting_area_square_metre", "optical_transmission", "observer_exposure_second")
    defs["photometry"] = operation("photometry.predict", {
        "source_model": {"const": "constant_rest_luminosity_rectangular_band"},
        "propagation": {"const": "isotropic_luminosity_distance_standard_redshift"},
        "constants_id": {"const": "si_2019_radiometric_definitions"},
        "inputs": array(obj({field: number for field in photometry_fields}), 65536),
        "requested_outputs": outputs(abi["photometry_tags"]["output"]),
        "resource_policy": obj({"maximum_rows": {**uint, "maximum": 65536},
            "maximum_native_bytes": {**uint, "maximum": 1073741824}})})
    bounded_text = {"type": "string", "minLength": 1, "maxLength": 1024}
    pool_meta = {"id": {"type": "string", "minLength": 1, "maxLength": 256},
        "source_role": enum(("measured", "fitted_summary", "calibration_asset", "synthetic_control")),
        "provenance": bounded_text}
    defs["sampled_photometry"] = operation("photometry.predict", {
        "source_model": {"const": "piecewise_linear_rest_luminosity_observed_optical_passband"},
        "propagation": {"const": "isotropic_luminosity_distance_standard_redshift"},
        "constants_id": {"const": "si_2019_radiometric_definitions"},
        "spectra": array(obj({**pool_meta, "rest_wavelength_metre": array(number,65536),
            "luminosity_watt_per_metre": array(number,65536)}),4096),
        "passbands": array(obj({**pool_meta, "calibration": enum(("fixed", "declared_uncertainty_excluded")),
            "calibration_provenance": bounded_text, "observed_wavelength_metre": array(number,65536),
            "optical_transmission": array(number,65536)}),4096),
        "exposures": array(obj({"spectrum_index": uint, "passband_index": uint,
            **{field:number for field in ("luminosity_distance_metre","redshift",
                "collecting_area_square_metre","observer_exposure_second")}}),65536),
        "requested_outputs": {**outputs(abi["photometry_tags"]["output"]), "maxItems":3},
        "resource_policy": obj({**{field:{**uint,"maximum":65536} for field in (
            "maximum_rows","maximum_total_samples","maximum_samples","maximum_segments")},
            "maximum_native_bytes": {**uint,"maximum":1073741824}})})
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
                                          "background", "statistics", "gaussian_posterior", "gaussian_predictive", "supernova", "bao", "bao_thermal", "sound_horizon", "photometry", "sampled_photometry")],
        "description": "Structural compiled requests; native domains, numerical gates and scientific qualifications remain distinct."}
