from flask import Flask, render_template
import requests
from datetime import datetime,  timedelta
from apscheduler.schedulers.background import BackgroundScheduler


def metro_times():
    # Define the GraphQL query
    stops={"Blindern":6332,"Kalbakken":5810}
    query = """
    {
    stopPlace(id: "NSR:StopPlace:##STOP_ID##") {
        name
        estimatedCalls(timeRange: 7200, numberOfDepartures: 10) {
        realtime
        aimedDepartureTime
        expectedDepartureTime
        destinationDisplay {
            frontText
        }
        serviceJourney {
            line {
            publicCode
            transportMode
            }
        }
        }
    }
    }
    """

    # Set the endpoint and headers
    url = "https://api.entur.io/journey-planner/v3/graphql"
    headers = {
        "Content-Type": "application/json",
        "ET-Client-Name": "your-client-id"  # Use your own identifier
    }

    # Send the request
    responses = {}
    departures = {}
    for i in stops:
        responses[i] = requests.post(url, json={"query": query.replace('##STOP_ID##',str(stops[i]))}, headers=headers)
        departures[i]=[]
    # Print result
    for i,name in enumerate(responses):
        if responses[name].status_code == 200:
            data = responses[name].json()
            
            for call in data["data"]["stopPlace"]["estimatedCalls"]:
                date_format = '%Y-%m-%dT%H:%M:%S%z'
                time=datetime.strptime(call["aimedDepartureTime"],date_format)
                time_expected=datetime.strptime(call["expectedDepartureTime"],date_format)
                delay=time_expected-time
                #departures[list(stops.keys())[i]].append(f"Line {call['serviceJourney']['line']['publicCode']:<4}  {call['destinationDisplay']['frontText']:<30} at {time.strftime('%H:%M')} expected at {time_expected.strftime('%H:%M')} delay {int(delay.total_seconds()//60):>3} minutes")
                departures[list(stops.keys())[i]].append({"Line":f"Line {call['serviceJourney']['line']['publicCode']}", 
                                                          "Destination":call['destinationDisplay']['frontText'],
                                                          "Time":time.strftime('%H:%M'),
                                                          "Expected":time_expected.strftime('%H:%M'),
                                                          "Delay":f"{int(delay.total_seconds()//60)} minutes"})
                departures[list(stops.keys())[i]]=departures[list(stops.keys())[i]][:4]
            
            #print(f"Line {call['serviceJourney']['line']['publicCode']} to {call['destinationDisplay']['frontText']} at {call['expectedDepartureTime']}")
            
            #print(f"Line {call['serviceJourney']['line']['publicCode']} to {call['destinationDisplay']['frontText']} at {call['aimedDepartureTime']} expected at  {call['expectedDepartureTime']}")
        else:
            print("FError:", responses[i].status_code)
    return(responses,str(delay),departures,stops)







app = Flask(__name__)

@app.route('/')
def index():
    responses,delay,departures,stops=metro_times()
    return render_template('./index.html', data1=departures['Blindern'], data2=departures['Kalbakken'])



if __name__ == '__main__':
    

    app.run(host='0.0.0.0', port=5000)
    
    
