#!/bin/bash
gunicorn -c gunicorn_config.py thmanager.wsgi

