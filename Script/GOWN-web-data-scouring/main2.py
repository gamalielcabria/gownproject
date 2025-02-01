from selenium import webdriver
from selenium.webdriver.chrome.service import Service

service = Service(executable_path='/home/glbcabria/Packages/installers/chromedriver-linux64/chromedriver')
options = webdriver.ChromeOptions()
driver = webdriver.Chromeget('https://quotes.toscrape.com/js')

element = driver.find_element_by_class_name("author")
# element = driver.find_element_by_tag_name("small")
# element = driver.find_element_by_xpath("//abc")

print(element.text)

driver.quit