from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.support.wait import WebDriverWait  

browser = webdriver.Edge()
try :
    browser.get ('https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRX5443883&o=acc_s%3Aa')

    print(browser.page_source)   #网页源代码
finally:
    browser.close()