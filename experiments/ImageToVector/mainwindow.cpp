#include "mainwindow.h"
#include "ui_mainwindow.h"

void msgbox(QString msg="message"){QMessageBox::information(0,"Message", msg);}

#include <QFileDialog>
#include <QPainter>

#include "imagetovector.h"
#include "PathFitter.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    connect(ui->inputImage, &QLabel::linkActivated, [&]{
        auto filenames = QFileDialog::getOpenFileNames(0,"Open image","","*.png");
        if(filenames.empty()) return;
        auto img = QImage(filenames.at(0));
        img.convertToFormat(QImage::Format_ARGB32_Premultiplied);
        ui->inputImage->setPixmap(QPixmap::fromImage(img));
    });

    connect(ui->convertButton, &QPushButton::clicked, [&]{
        // Try to fill either side of a boundary and select one that doesn't make an outside mess
        auto grid = make_grid<double>( ui->inputImage->pixmap()->toImage() );
        auto start = first_nonempty(grid);
        auto neighbours = empty_nighbours(grid, start);
        std::array<size_t,2> selected = {-1,-1};
        for(auto cell : neighbours)
        {
            auto copy_grid = grid;
            floodfill_stacked(cell[0], cell[1], copy_grid);

            if(copy_grid[0][0] == 0.0) continue;
            selected = cell;
            break;
        }
        floodfill_stacked(selected[0], selected[1], grid);

        // Generate filled image
        auto filled_img = ui->inputImage->pixmap()->toImage();
        for(int x = 0; x < filled_img.width(); x++){
            for(int y = 0; y < filled_img.height(); y++){
                auto v = grid[y][x] * 255;
                filled_img.setPixel(x, y, qRgb(v,v,v));
            }
        }

        // Collect any possible two boundaries:
        typedef std::vector< std::vector<double> > GridType;
        QVector<  std::vector< std::array<double,2> > > boundaries;
        boundaries << MarchingSquares< std::array<double,2>, GridType >::march(grid, 0);
		floodfill_stacked(0, 0, grid);
        boundaries << MarchingSquares< std::array<double,2>, GridType >::march(grid, 1);

		QPainterPath letterPath;

        // Simplify each boundary by fitting
        for(auto b : boundaries)
        {
			if (b.empty()) continue;

            if(ui->vectorOption->isChecked())
            {
                std::vector<Vec2f> pnts;
                for(auto coord : b) pnts.push_back(Vec2f(coord[0], coord[1]));

                double error = 2;
                auto path = PathFitter::FitCurve(pnts, 0, pnts.size(), error);
                letterPath.addPath(path);
            }

            if(ui->smoothingOption->isChecked())
            {
                QPolygonF poly;
                for(auto coord : b) poly.push_back(QPoint(coord[0], coord[1]));

                poly = resamplePolygon(poly, b.size() * 0.9);
                poly = smoothPolygon(poly, 3);

                letterPath.addPolygon(poly);
            }

            if(ui->nopostOption->isChecked())
            {
                QPolygonF poly;
                for(auto coord : b) poly.push_back(QPoint(coord[0], coord[1]));
                letterPath.addPolygon(poly);
            }
        }

        // Visualize
        QPainter painter(&filled_img);
        painter.setPen(QPen(Qt::red,2));
        painter.drawPath(letterPath);

        // Update
        ui->inputImage->setPixmap(QPixmap::fromImage(filled_img));

        // Output vectorized version
        auto size = ui->outputImage->size();
        QImage vector_img(size.width(), size.height(), QImage::Format_ARGB32_Premultiplied);
        vector_img.fill(Qt::white);
        QPainter vec_painter(&vector_img);
        vec_painter.translate(vector_img.rect().center());
        auto letterRect = letterPath.boundingRect();
        letterPath.translate(-letterRect.center());

        double s = std::min(double(size.height()) / double(letterRect.height()) * 0.8,
                            double(size.width()) / double(letterRect.width()) * 0.8);

        vec_painter.scale( s, s );
        vec_painter.fillPath(letterPath, Qt::black);
        ui->outputImage->setPixmap(QPixmap::fromImage(vector_img));
    });

    connect(ui->saveButton, &QPushButton::clicked, [&]{
        int counter = 0;
        QString filename = ui->lineEdit->text() + QString::number(counter) + ".png";
        ui->outputImage->pixmap()->toImage().save(QFileDialog::getSaveFileName(0,"Open image",filename,"*.png"));
    });

	// This is just repeating above codes.. todo: turn to function
	connect(ui->batchButton, &QPushButton::clicked, [&]{
		auto filenames = QFileDialog::getOpenFileNames(0, "Open image", "", "*.png");
		if (filenames.empty()) return;

		for (auto f : filenames)
		{
			auto img = QImage(f);
			img.convertToFormat(QImage::Format_ARGB32_Premultiplied);
			ui->inputImage->setPixmap(QPixmap::fromImage(img));

			// Try to fill either side of a boundary and select one that doesn't make an outside mess
			auto grid = make_grid<double>(ui->inputImage->pixmap()->toImage());
			auto start = first_nonempty(grid);
			auto neighbours = empty_nighbours(grid, start);
			std::array<size_t, 2> selected = { -1, -1 };
			for (auto cell : neighbours)
			{
				auto copy_grid = grid;
				floodfill_stacked(cell[0], cell[1], copy_grid);

				if (copy_grid[0][0] == 0.0) continue;
				selected = cell;
				break;
			}
			floodfill_stacked(selected[0], selected[1], grid);

			// Generate filled image
			auto filled_img = ui->inputImage->pixmap()->toImage();
			for (int x = 0; x < filled_img.width(); x++){
				for (int y = 0; y < filled_img.height(); y++){
					auto v = grid[y][x] * 255;
					filled_img.setPixel(x, y, qRgb(v, v, v));
				}
			}

			// Collect any possible two boundaries:
			typedef std::vector< std::vector<double> > GridType;
			QVector<  std::vector< std::array<double, 2> > > boundaries;
			boundaries << MarchingSquares< std::array<double, 2>, GridType >::march(grid, 0);
			floodfill_stacked(0, 0, grid);
			boundaries << MarchingSquares< std::array<double, 2>, GridType >::march(grid, 1);

			QPainterPath letterPath;

			// Simplify each boundary by fitting
			for (auto b : boundaries)
			{
				if (b.empty()) continue;

				if (ui->vectorOption->isChecked())
				{
					std::vector<Vec2f> pnts;
					for (auto coord : b) pnts.push_back(Vec2f(coord[0], coord[1]));

					double error = 2;
					auto path = PathFitter::FitCurve(pnts, 0, pnts.size(), error);
					letterPath.addPath(path);
				}

				if (ui->smoothingOption->isChecked())
				{
					QPolygonF poly;
					for (auto coord : b) poly.push_back(QPoint(coord[0], coord[1]));

					poly = resamplePolygon(poly, b.size() * 0.9);
					poly = smoothPolygon(poly, 3);

					letterPath.addPolygon(poly);
				}

				if (ui->nopostOption->isChecked())
				{
					QPolygonF poly;
					for (auto coord : b) poly.push_back(QPoint(coord[0], coord[1]));
					letterPath.addPolygon(poly);
				}
			}

			// Visualize
			QPainter painter(&filled_img);
			painter.setPen(QPen(Qt::red, 2));
			painter.drawPath(letterPath);

			// Update
			ui->inputImage->setPixmap(QPixmap::fromImage(filled_img));

			// Output vectorized version
			auto size = ui->outputImage->size();
			QImage vector_img(size.width(), size.height(), QImage::Format_ARGB32_Premultiplied);
			vector_img.fill(Qt::white);
			QPainter vec_painter(&vector_img);
			vec_painter.translate(vector_img.rect().center());
			auto letterRect = letterPath.boundingRect();
			letterPath.translate(-letterRect.center());

			double s = std::min(double(size.height()) / double(letterRect.height()) * 0.8,
				double(size.width()) / double(letterRect.width()) * 0.8);

			vec_painter.scale(s, s);
			vec_painter.fillPath(letterPath, Qt::black);
			ui->outputImage->setPixmap(QPixmap::fromImage(vector_img));

			// Save output:
			vector_img.save(QString("%1.converted.png").arg(f));
		}
	});
}

MainWindow::~MainWindow()
{
    delete ui;
}
