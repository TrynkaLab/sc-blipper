#!/usr/bin/env nextflow
//include { validateParameters } from 'plugin/nf-schema'

// Validate parameters plugin
//validateParameters()

// Run pipeline
include { cnmf } from './workflows/run_cnmf.nf'