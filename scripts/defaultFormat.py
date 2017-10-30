#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import re, datetime

from googleapiclient.discovery import build
from oauth2client.service_account import ServiceAccountCredentials

# connect to api
api = build('sheets', 'v4',
    credentials=ServiceAccountCredentials.from_json_keyfile_name(
        'client_secret.json', 'https://www.googleapis.com/auth/spreadsheets')
    ).spreadsheets()

# sheet ID
ID = '1NEta9svTbKyTx48UdQIyLW2J6TWGrNwqvT6m1mjqQFQ'

api.batchUpdate(spreadsheetId=ID, body = { 'requests': [
  { 'updateSpreadsheetProperties': {
    'properties': {
      'defaultFormat': {
        'padding': { 'top': 0, 'right': 1, 'bottom': 0, 'left': 1 }
      },
    }, 'fields': 'defaultFormat(padding)'
  } }
]}).execute()
