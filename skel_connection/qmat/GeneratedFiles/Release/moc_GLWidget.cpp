/****************************************************************************
** Meta object code from reading C++ file 'GLWidget.h'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.9.9)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../../GLWidget.h"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'GLWidget.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.9.9. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
QT_WARNING_PUSH
QT_WARNING_DISABLE_DEPRECATED
struct qt_meta_stringdata_GLWidget_t {
    QByteArrayData data[22];
    char stringdata0[251];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_GLWidget_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_GLWidget_t qt_meta_stringdata_GLWidget = {
    {
QT_MOC_LITERAL(0, 0, 8), // "GLWidget"
QT_MOC_LITERAL(1, 9, 16), // "xRotationChanged"
QT_MOC_LITERAL(2, 26, 0), // ""
QT_MOC_LITERAL(3, 27, 5), // "angle"
QT_MOC_LITERAL(4, 33, 16), // "yRotationChanged"
QT_MOC_LITERAL(5, 50, 16), // "zRotationChanged"
QT_MOC_LITERAL(6, 67, 11), // "ZoomChanged"
QT_MOC_LITERAL(7, 79, 6), // "factor"
QT_MOC_LITERAL(8, 86, 16), // "TranslateChanged"
QT_MOC_LITERAL(9, 103, 2), // "dx"
QT_MOC_LITERAL(10, 106, 2), // "dy"
QT_MOC_LITERAL(11, 109, 20), // "setXYRotationChanged"
QT_MOC_LITERAL(12, 130, 6), // "xAngle"
QT_MOC_LITERAL(13, 137, 6), // "yAngle"
QT_MOC_LITERAL(14, 144, 20), // "setXZRotationChanged"
QT_MOC_LITERAL(15, 165, 6), // "zAngle"
QT_MOC_LITERAL(16, 172, 18), // "setRotationChanged"
QT_MOC_LITERAL(17, 191, 12), // "setXRotation"
QT_MOC_LITERAL(18, 204, 12), // "setYRotation"
QT_MOC_LITERAL(19, 217, 12), // "setZRotation"
QT_MOC_LITERAL(20, 230, 7), // "setZoom"
QT_MOC_LITERAL(21, 238, 12) // "setTranslate"

    },
    "GLWidget\0xRotationChanged\0\0angle\0"
    "yRotationChanged\0zRotationChanged\0"
    "ZoomChanged\0factor\0TranslateChanged\0"
    "dx\0dy\0setXYRotationChanged\0xAngle\0"
    "yAngle\0setXZRotationChanged\0zAngle\0"
    "setRotationChanged\0setXRotation\0"
    "setYRotation\0setZRotation\0setZoom\0"
    "setTranslate"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_GLWidget[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
      13,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       5,       // signalCount

 // signals: name, argc, parameters, tag, flags
       1,    1,   79,    2, 0x06 /* Public */,
       4,    1,   82,    2, 0x06 /* Public */,
       5,    1,   85,    2, 0x06 /* Public */,
       6,    1,   88,    2, 0x06 /* Public */,
       8,    2,   91,    2, 0x06 /* Public */,

 // slots: name, argc, parameters, tag, flags
      11,    2,   96,    2, 0x0a /* Public */,
      14,    2,  101,    2, 0x0a /* Public */,
      16,    3,  106,    2, 0x0a /* Public */,
      17,    1,  113,    2, 0x0a /* Public */,
      18,    1,  116,    2, 0x0a /* Public */,
      19,    1,  119,    2, 0x0a /* Public */,
      20,    1,  122,    2, 0x0a /* Public */,
      21,    2,  125,    2, 0x0a /* Public */,

 // signals: parameters
    QMetaType::Void, QMetaType::Int,    3,
    QMetaType::Void, QMetaType::Int,    3,
    QMetaType::Void, QMetaType::Int,    3,
    QMetaType::Void, QMetaType::Double,    7,
    QMetaType::Void, QMetaType::Double, QMetaType::Double,    9,   10,

 // slots: parameters
    QMetaType::Void, QMetaType::Int, QMetaType::Int,   12,   13,
    QMetaType::Void, QMetaType::Int, QMetaType::Int,   12,   15,
    QMetaType::Void, QMetaType::Int, QMetaType::Int, QMetaType::Int,   12,   13,   15,
    QMetaType::Void, QMetaType::Int,    3,
    QMetaType::Void, QMetaType::Int,    3,
    QMetaType::Void, QMetaType::Int,    3,
    QMetaType::Void, QMetaType::Double,    7,
    QMetaType::Void, QMetaType::Double, QMetaType::Double,    9,   10,

       0        // eod
};

void GLWidget::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        GLWidget *_t = static_cast<GLWidget *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->xRotationChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 1: _t->yRotationChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 2: _t->zRotationChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 3: _t->ZoomChanged((*reinterpret_cast< double(*)>(_a[1]))); break;
        case 4: _t->TranslateChanged((*reinterpret_cast< double(*)>(_a[1])),(*reinterpret_cast< double(*)>(_a[2]))); break;
        case 5: _t->setXYRotationChanged((*reinterpret_cast< int(*)>(_a[1])),(*reinterpret_cast< int(*)>(_a[2]))); break;
        case 6: _t->setXZRotationChanged((*reinterpret_cast< int(*)>(_a[1])),(*reinterpret_cast< int(*)>(_a[2]))); break;
        case 7: _t->setRotationChanged((*reinterpret_cast< int(*)>(_a[1])),(*reinterpret_cast< int(*)>(_a[2])),(*reinterpret_cast< int(*)>(_a[3]))); break;
        case 8: _t->setXRotation((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 9: _t->setYRotation((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 10: _t->setZRotation((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 11: _t->setZoom((*reinterpret_cast< double(*)>(_a[1]))); break;
        case 12: _t->setTranslate((*reinterpret_cast< double(*)>(_a[1])),(*reinterpret_cast< double(*)>(_a[2]))); break;
        default: ;
        }
    } else if (_c == QMetaObject::IndexOfMethod) {
        int *result = reinterpret_cast<int *>(_a[0]);
        {
            typedef void (GLWidget::*_t)(int );
            if (*reinterpret_cast<_t *>(_a[1]) == static_cast<_t>(&GLWidget::xRotationChanged)) {
                *result = 0;
                return;
            }
        }
        {
            typedef void (GLWidget::*_t)(int );
            if (*reinterpret_cast<_t *>(_a[1]) == static_cast<_t>(&GLWidget::yRotationChanged)) {
                *result = 1;
                return;
            }
        }
        {
            typedef void (GLWidget::*_t)(int );
            if (*reinterpret_cast<_t *>(_a[1]) == static_cast<_t>(&GLWidget::zRotationChanged)) {
                *result = 2;
                return;
            }
        }
        {
            typedef void (GLWidget::*_t)(double );
            if (*reinterpret_cast<_t *>(_a[1]) == static_cast<_t>(&GLWidget::ZoomChanged)) {
                *result = 3;
                return;
            }
        }
        {
            typedef void (GLWidget::*_t)(double , double );
            if (*reinterpret_cast<_t *>(_a[1]) == static_cast<_t>(&GLWidget::TranslateChanged)) {
                *result = 4;
                return;
            }
        }
    }
}

const QMetaObject GLWidget::staticMetaObject = {
    { &QGLWidget::staticMetaObject, qt_meta_stringdata_GLWidget.data,
      qt_meta_data_GLWidget,  qt_static_metacall, nullptr, nullptr}
};


const QMetaObject *GLWidget::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *GLWidget::qt_metacast(const char *_clname)
{
    if (!_clname) return nullptr;
    if (!strcmp(_clname, qt_meta_stringdata_GLWidget.stringdata0))
        return static_cast<void*>(this);
    return QGLWidget::qt_metacast(_clname);
}

int GLWidget::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QGLWidget::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 13)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 13;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 13)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 13;
    }
    return _id;
}

// SIGNAL 0
void GLWidget::xRotationChanged(int _t1)
{
    void *_a[] = { nullptr, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}

// SIGNAL 1
void GLWidget::yRotationChanged(int _t1)
{
    void *_a[] = { nullptr, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 1, _a);
}

// SIGNAL 2
void GLWidget::zRotationChanged(int _t1)
{
    void *_a[] = { nullptr, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 2, _a);
}

// SIGNAL 3
void GLWidget::ZoomChanged(double _t1)
{
    void *_a[] = { nullptr, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 3, _a);
}

// SIGNAL 4
void GLWidget::TranslateChanged(double _t1, double _t2)
{
    void *_a[] = { nullptr, const_cast<void*>(reinterpret_cast<const void*>(&_t1)), const_cast<void*>(reinterpret_cast<const void*>(&_t2)) };
    QMetaObject::activate(this, &staticMetaObject, 4, _a);
}
QT_WARNING_POP
QT_END_MOC_NAMESPACE
