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

info = api.get(spreadsheetId=ID).execute()
sheets = [ s['properties']['title'] for s in info['sheets'] ]

today = datetime.date.today().strftime('%d %b %Y')
title = today

# add suffix to new title if needed
def suffix(s):
  if not s.startswith(title): return None
  s = s[len(title):]
  if len(s)==0: return 0
  try: return int(re.split('\((\d+)\)$',s)[1])
  except: return None

dups = filter(lambda x: x is not None, [ suffix(s) for s in sheets ])
if len(dups)!=0:
  title += ' ({})'.format(max(dups)+1)

# add sheet
# https://developers.google.com/sheets/api/samples/sheet#add_a_sheet
print "Adding sheet: ", title
result = api.batchUpdate(spreadsheetId=ID, body = \
{ 'requests': [
  { 'addSheet': {
      'properties': {
        'title': title,
        'gridProperties': { 'rowCount': 1000, 'columnCount': 11 }
      }
    }
  }]
}).execute()

print result
sheetId = result['replies'][0]['addSheet']['properties']['sheetId']

# header
api.values().update(spreadsheetId=ID, range=title+'!A1', body = \
{ 'values': [
['','','','','','unc','','unc','signif','','reco'],
['Variable','[',')','width','sig','√(Σ(s²))','bkg','√(Σ(b²))','s/√(s+b)','s/(s+b)','purity']
]}, valueInputOption='RAW').execute()

# conditional formatting
def cfmt(rule_t,ranges,rule,i=0):
    return { 'addConditionalFormatRule': {
        'rule': { 'ranges': ranges, rule_t : rule }, 'index': i
      } }

ranges = [
  { 'sheetId': sheetId,
    'startRowIndex': 2, 'startColumnIndex': 10, 'endColumnIndex': 11 }
]
api.batchUpdate(spreadsheetId=ID, body = \
{ 'requests': [
  cfmt('booleanRule',ranges,{
    'condition': {
      'type': 'NUMBER_LESS', 'values': [ { 'userEnteredValue': '0.4' } ]
    },
    'format': {
      'textFormat': {
        'bold': True, 'foregroundColor': { 'red': 204./255 }
      }
    }
  },0),
  cfmt('booleanRule',ranges,{
    'condition': {
      'type': 'NUMBER_LESS', 'values': [ { 'userEnteredValue': '0.5' } ]
    },
    'format': {
      'textFormat': {
        'bold': True, 'foregroundColor': { 'red': 255./255, 'green': 102./255 }
      }
    }
  },1),
  cfmt('booleanRule',ranges,{
    'condition': {
      'type': 'NUMBER_LESS', 'values': [ { 'userEnteredValue': '0.75' } ]
    },
    'format': {
      'textFormat': {
        'bold': True, 'foregroundColor': { 'blue': 153./255 }
      }
    }
  },2),
  cfmt('booleanRule',ranges,{
    'condition': {
      'type': 'NUMBER_GREATER_THAN_EQ', 'values': [ { 'userEnteredValue': '0.75' } ]
    },
    'format': {
      'textFormat': {
        'bold': True, 'foregroundColor': { 'green': 102./255 }
      }
    }
  },3)
]}).execute()
# 'gradientRule': {
#   'minpoint': { 'color': { 'red': float.fromhex('0x0.CC') },
#                 'type': 'NUMBER', 'value': '0' },
#   'midpoint': { 'color': { 'blue': float.fromhex('0x0.99') },
#                 'type': 'NUMBER', 'value': '0.7' },
#   'maxpoint': { 'color': { 'green': float.fromhex('0x0.66') },
#                 'type': 'NUMBER', 'value': '1' }
# }

