"""Routine engineering coverage; every other registered test is opt-in regression.

Membership is reviewed by purpose, not by a timing threshold. These tests exercise
ABI descriptors, ownership, arithmetic admission, small algebra, and refusals.
Peer/reference, cosmological solvers, recovery and scientific consumers remain full.
"""
import re

NATIVE_FAST = {
    **{f'{name}_contract': f'test_{name}' for name in (
        'abi', 'abi_hostile', 'quantities', 'quantities_hostile',
        'quantity_abi_hostile', 'numerics', 'arithmetic_profile', 'numerics_hostile',
        'observations', 'numerical_abi_hostile', 'observation_abi_hostile',
        'observations_hostile', 'statistics', 'statistics_hostile',
        'gaussian_abi_hostile', 'profile_operator_hostile',
        'observation_ownership', 'observation_ownership_hostile')},
    'fits_codec_disabled_contract': 'test_fits_reader_disabled',
    'lifetime_scale_hostile': 'test_lifetime_scale_hostile',
}
# Immutable records/fault transport, decimal JSON, quantities, scalar algebra,
# synthetic ingestion, operation-owned outcomes and duplicate-key/session refusal.
CLI_FAST = ('cli', 'json_precision', 'outcome_cli', 'quantities_cli',
            'numerics_cli', 'ingestion', 'strict_json_cli')


def native_build_targets(suite):
    if suite == 'core':
        return ['irred_core']
    if suite == 'fast':
        return sorted(set(NATIVE_FAST.values()))
    if suite == 'full':
        return []  # CMake's complete default target; no test is deleted.
    raise ValueError(f'unknown native suite: {suite}')


def classify_native(inventory):
    names = [t['name'] for t in inventory['tests']]
    if len(names) != len(set(names)):
        raise ValueError('duplicate registered native test')
    missing = set(NATIVE_FAST) - set(names)
    if missing:
        raise ValueError(f'routine native tests missing: {sorted(missing)}')
    return {'fast': sorted(NATIVE_FAST), 'full_only': sorted(set(names) - set(NATIVE_FAST))}


def native_regex():
    return '^(' + '|'.join(re.escape(n) for n in sorted(NATIVE_FAST)) + ')$'


def cli_targets(available, suite):
    available = set(available)
    missing = set(CLI_FAST) - available
    if missing:
        raise ValueError(f'routine CLI targets missing: {sorted(missing)}')
    return sorted(available if suite == 'full' else CLI_FAST)


def installed_sources(root, suite):
    if suite == 'fast':
        return [root / 'cpp/tests/test_abi.cpp']
    if suite == 'full':
        return [root / name for name in (
            'cpp/tests/test_installed_consumer.cpp',
            'cpp/tests/test_installed_continuous_cmb_projection.cpp',
            'cpp/examples/windowed_linear_power.cpp',
            'cpp/tests/test_effective_fluid_installed.cpp',
            'cpp/tests/test_two_deflector_forward_installed.cpp',
            'cpp/tests/test_installed_finite_opacity_source.cpp',
            'cpp/tests/test_installed_ideal_acoustic.cpp',
            'cpp/examples/quintessence.cpp',
            'cpp/tests/test_dgp_growth_installed.cpp',
            'cpp/tests/test_nfw_halo_installed.cpp',
            'cpp/examples/decaying_matter.cpp',
            'cpp/tests/test_curved_flrw_sdk.cpp')]
    raise ValueError(f'unknown installed suite: {suite}')
