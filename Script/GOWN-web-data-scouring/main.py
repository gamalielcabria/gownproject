import requests
from  requests_html import HTMLSession
from bs4 import BeautifulSoup

url = 'https://groundwatermonitoring.alberta.ca/#/overview/activewells/station/41425/Vermilion%202597E_0294/info'

session = HTMLSession()
response = session.get(url)
response.html.render()

print(response.html)
print(response.html.html)

# soup = BeautifulSoup(response.content, 'html.parser')

# print(soup.find_all('wwp-station-info-albg', class_='lg-screen LANDSCAPE'))