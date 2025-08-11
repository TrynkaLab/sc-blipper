#!/usr/bin/env nextflow
//include { validateParameters } from 'plugin/nf-schema'

// Validate parameters plugin
//validateParameters()

// Run pipeline
include { run_cnmf } from './workflows/run_cnmf.nf'