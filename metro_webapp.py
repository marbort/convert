from flask import Flask

app = Flask(__name__)
@app.route('/')
def hello():
    return "Hello, World!"
@app.route('/goodbye')
def goodbye():
    return "Goodbye, World!"
@app.route('/greet/<name>')
def greet(name):
    return f"Hello, {name}!"


if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5000)
# To run the application, save this code in a file named app.py and run it using the command:
# python app.py
# The application will be accessible at http://localhost:5000
# You can test the different routes by visiting:
# - http://localhost:5000/ for "Hello, World!"
# - http://localhost:5000/goodbye for "Goodbye, World!"
# - http://localhost:5000/greet/YourName for "Hello, YourName!"
# Make sure to replace YourName with any name you want to greet.