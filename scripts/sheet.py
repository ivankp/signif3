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
# print info

title = datetime.date.today().strftime('%d %b %Y')

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
result = api.batchUpdate(spreadsheetId=ID, body = { 'requests': [
  { 'addSheet': {
      'properties': {
        'title': title,
        'gridProperties': { 'rowCount': 25, 'columnCount': 11 }
      }
  } }
]}).execute()

# print result
sheetId = result['replies'][0]['addSheet']['properties']['sheetId']
print 'sheetId:', sheetId

def col(i):
  return { 'sheetId': sheetId, 'startRowIndex': 2,
           'startColumnIndex': i, 'endColumnIndex': i+1 }

# conditional formatting
def cond_color(ranges,colors):
    return [ {
      'addConditionalFormatRule': {
        'rule': {
          'ranges': ranges,
          'booleanRule': {
            'condition': {
              'type': 'NUMBER_LESS',
              'values': [ { 'userEnteredValue': x[1][0] } ]
            },
            'format': { 'textFormat': { 'foregroundColor': x[1][1] } }
          }
        }, 'index': x[0]
      }
    } for x in enumerate(colors) ]

def color(c):
  c1 = divmod(c,256)
  c2 = divmod(c1[0],256)
  return { 'red': c2[0]/255., 'green': c2[1]/255., 'blue': c1[1]/255. }

def fmt(ranges,fields):
  return [ {
    'repeatCell': {
      'range': r,
      'cell': { 'userEnteredFormat': fields },
      'fields': 'userEnteredFormat(' + ','.join([ f for f in fields ]) +')'
    } } for r in ranges ]

def width(w,a,b):
  return { "updateDimensionProperties": {
    "range": {
      "sheetId": sheetId,
      "dimension": "COLUMNS",
      "startIndex": a,
      "endIndex": b
    }, "properties": { "pixelSize": w }, "fields": "pixelSize"
  } }

# https://developers.google.com/sheets/api/samples/formatting
api.batchUpdate(spreadsheetId=ID, body = { 'requests':
  [ { 'updateSheetProperties': {
        'properties': {
          'sheetId': sheetId,
          'gridProperties': { 'frozenRowCount': 2 }
        },
        'fields': 'gridProperties.frozenRowCount'
    } },
    width(52,1,11),
    { "updateDimensionProperties": {
        "range": { "sheetId": sheetId, "dimension": "ROWS" },
        "properties": { "pixelSize": 17 }, "fields": "pixelSize"
    } }
  ] + \
  fmt([ { 'sheetId': sheetId } ],
        { 'padding': { 'top': 0, 'right': 2, 'bottom': 0, 'left': 2 } }) + \
  fmt([ { 'sheetId': sheetId,
          'startRowIndex': 0, 'endRowIndex': 2,
          'startColumnIndex': 0, 'endColumnIndex': 11 } ],
        { 'horizontalAlignment' : 'CENTER' }) + \
  fmt([ col(3), col(5), col(7) ],
        { 'textFormat': { 'foregroundColor': color(0x434343) } }) + \
  fmt([ col(8), col(10) ],
        { 'textFormat': { 'bold': True } }) + \
  cond_color( col(8),
    [ ('1'  , color(0xCC0000)),
      ('2'  , color(0xFF6600)),
      ('2.3', color(0x000099)),
      ('100', color(0x006600))
    ]
  ) + \
  cond_color( col(10),
    [ ('0.4' , color(0xCC0000)),
      ('0.5' , color(0xFF6600)),
      ('0.75', color(0x000099)),
      ('1'   , color(0x006600))
    ]
  ) + \
  [ { 'addConditionalFormatRule': {
        'rule': {
          'ranges': col(0),
          'booleanRule': {
            'condition': {
              'type': 'CUSTOM_FORMULA',
              'values': [ { 'userEnteredValue':
                '=OR(ROW()=3,AND(ISBLANK(A2),ISBLANK(B2)))' } ]
            },
            'format': { 'textFormat': { 'bold': True } }
          }
        }, 'index': 0
      }
  } ]
}).execute()

# header
api.values().update(spreadsheetId=ID, range=title+'!A1', body = \
{ 'values': [
['','','','','','unc','','unc','signif','','reco'],
['Variable','[',')','width','sig','√(Σ s²)','bkg','√(Σ b²)','s/√(s+b)','s/(s+b)','purity']
]}, valueInputOption='RAW').execute()

