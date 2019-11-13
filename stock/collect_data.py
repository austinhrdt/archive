from iexfinance.stocks import get_historical_intraday
from pandas.tseries.offsets import BDay # B for buisiness, not birth
import pandas as pd
import tqdm as tqdm
import datetime
import time


# obtain intraday 
def get_intraday(ticker, date):
    df = pd.DataFrame(get_historical_intraday(ticker, date))
    if df.empty or df['marketAverage'].isnull().any().any() : return None
    df['date'] = pd.to_datetime(df['date'] + ' ' + df['minute'])
    return df.set_index('date')['marketAverage']

# obtain data for multiple days
def get_multiple_days(ticker, period):
    intraday_data = []
    for day in reversed(range(period)):
        date = datetime.datetime.today() - BDay(day)
        intraday_data.append(get_intraday(ticker, date))
    df = pd.concat(intraday_data).rename(ticker)
    return df[df != -1].dropna()

# combine data for each stock into a dataframe
def merge_data(period):
    stocks, tickers = [], get_tickers()
    for ticker in tickers:
        try: stocks.append(get_multiple_days(ticker, period))
        except: pass
    df = pd.concat(stocks, axis=1)
    return df

# returns a list of tickers.
def get_tickers():
    return list(pd.read_csv('tickers.csv')['Symbol'])



start = time

time.sleep(5)

end = time

print(end, start)
    
