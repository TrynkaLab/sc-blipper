#!/usr/bin/env nextflow
//include { validateParameters } from 'plugin/nf-schema'

// Validate parameters plugin
//validateParameters()

// Run pipeline
include { run_pipeline } from './workflows/run_pipeline.nf'